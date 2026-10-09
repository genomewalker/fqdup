#include "../src/merge_r8qkeep/atomic_family_union.hpp"
#include "../src/merge_r8qkeep/worker_pool.hpp"
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
using namespace fqdup::r8qkeep;
int main(){
    WorkerPool pool(8);std::mt19937 rng(2605);
    for(unsigned mode=0;mode<4;++mode){
        const int n=32000;std::vector<std::pair<int,int>> edges;
        for(int i=1;i<n;++i)if(i%17==0)edges.emplace_back(i,i-1);
        for(int i=0;i<120000;++i){int a=rng()%n,b=rng()%n;if(mode==0)b=a;else if(mode==1)b=std::min(n-1,a+1);else if(mode==2)b=a/128*128+b%128;edges.emplace_back(a,b);}
        std::vector<int> ref(n);std::iota(ref.begin(),ref.end(),0);
        const auto find=[&](int i){while(ref[i]!=i){ref[i]=ref[ref[i]];i=ref[i];}return i;};
        for(auto [a,b]:edges){a=find(a);b=find(b);ref[std::max(a,b)]=std::min(a,b);}
        for(int repeats=0;repeats<4;++repeats)for(int workers:{1,8}){
            std::vector<int> parent(n);std::iota(parent.begin(),parent.end(),0);
            for(int i=1;i<n;++i)if(i%17==0)parent[i]=i-1;
            AtomicFamilyUnion sets(parent);std::vector<std::future<void>> jobs;
            for(int w=0;w<workers;++w)jobs.push_back(pool.submit([&,w]{for(std::size_t i=w;i<edges.size();i+=workers){const auto e=edges[i];sets.unite(e.first,e.second);}}));
            for(auto& job:jobs)job.get();
            for(int i=0;i<n;++i)if(sets.root(i)!=find(i)||parent[i]>i)throw std::runtime_error("atomic component parity");
        }
    }
    std::cout<<"atomic monotone family union PASS\n";
}
