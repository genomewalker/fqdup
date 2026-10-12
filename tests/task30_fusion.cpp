#include "../src/merge_r8qkeep/repair.cpp"
#define fqalign task30_native_alignment
#include "../src/merge_r8qkeep/native_headers/fqdup/pair_quality/native_quality.hpp"
#undef fqalign
#include <random>
#include <future>

namespace repair = fqdup_r8qkeep_repair_internal::fqalign;
namespace native = task30_native_alignment;
namespace cache = fqdup_shared_alignment;
static void require(bool value) {
    if (!value) throw std::runtime_error("Task30 cross-stage alignment mismatch");
}
static void same(const cache::Correspondence& a, const cache::Correspondence& b) {
    require(a.best==b.best && a.paths==b.paths && a.eligible==b.eligible &&
        a.reason==b.reason && a.invariant==b.invariant && a.certified==b.certified);
}
static void same(const cache::Geometry& a, const cache::Geometry& b) {
    require(a.insert==b.insert && a.informative==b.informative &&
        a.mismatches==b.mismatches && a.fraction==b.fraction &&
        a.second==b.second && a.eligible==b.eligible);
}
static void exercise(unsigned seed) {
    std::mt19937 rng(seed);
    native::BandedCorrespondence original;
    repair::BandedCorrespondence optimized;
    for (int length : {0,1,11,20,31,60,100,150,511,512,513}) {
        for (int variant=0; variant<6; ++variant) {
            std::string a(length,'A');
            for (char& c:a) c="ACGT"[rng()%4];
            auto b=a;
            if (variant==1) for (char& c:b) c="ACGT"[rng()%4];
            if (variant==2 && length>7) b[7]='N';
            if (variant==3 && length>20) b.insert(b.begin()+length/2,'G');
            if (variant==4) { a.assign(length,'N'); b=a; }
            if (variant==5 && length>4) b.erase(b.begin()+length/2);
            for (int insert : {-1,0,length/2,length,length+7}) {
                const auto expected=original.align_uncached(a,b,insert);
                same(expected,optimized.align_uncached(a,b,insert));
                cache::Pair products;
                cache::Scope scope(&products);
                same(expected,original.align(a,b,insert));
                require(products.misses==1);
                same(expected,optimized.align(a,b,insert));
                require(products.misses==1); // Native product, not a second DP.
                auto copy=optimized.align(a,b,insert); copy.certified.clear();
                same(expected,original.align(a,b,insert)); // Immutable entry.
                optimized.align(a,b,insert+1); require(products.misses==2);
                optimized.align(a+"A",b,insert); require(products.misses==3);
                cache::Pair other;
                { cache::Scope nested(&other); optimized.align(a,b,insert); }
                require(other.misses==1);
                same(expected,optimized.align(a,b,insert));
                require(products.misses==3);
            }
            const auto raw_b=native::reverse_complement(b);
            const auto expected=native::quality_independent_geometry_uncached(a,raw_b);
            same(expected,repair::quality_independent_geometry_uncached(a,raw_b));
            cache::Pair products; cache::Scope scope(&products);
            same(expected,native::quality_independent_geometry(a,raw_b));
            same(expected,repair::quality_independent_geometry(a,raw_b));
            require(products.misses==1);
        }
    }
    require(cache::active==nullptr);
}
int main() {
    exercise(3004);
    auto worker=std::async(std::launch::async,[]{exercise(3005);});
    exercise(3006); worker.get();
    std::cout << "TASK30_NATIVE_REPAIR_PRODUCTS_PASS\n";
}
