#pragma once
#include <algorithm>
#include <vector>
namespace fqdup::r8qkeep {
// During the parallel phase every parent access uses an atomic builtin. Each
// successful link lowers its root ID, so links cannot form a cycle. A failed
// CAS retries both roots; no accepted graph edge is lost to a competing union.
// Initial exact-duplicate links must also point to smaller IDs. Once workers
// join, ordinary serial path compression may resume on this same flat array.
class AtomicFamilyUnion {
    std::vector<int>& parent_;
public:
    explicit AtomicFamilyUnion(std::vector<int>& parent):parent_(parent){}
    int root(int i)const {
        for(;;){const int p=__atomic_load_n(&parent_[i],__ATOMIC_RELAXED);if(p==i)return i;
            const int next=__atomic_load_n(&parent_[p],__ATOMIC_RELAXED);
            if(next!=p){int expected=p;__atomic_compare_exchange_n(&parent_[i],&expected,next,false,__ATOMIC_RELAXED,__ATOMIC_RELAXED);}
            i=p;
        }
    }
    void unite(int a,int b)const {
        for(;;){a=root(a);b=root(b);if(a==b)return;if(a<b)std::swap(a,b);
            int expected=a;if(__atomic_compare_exchange_n(&parent_[a],&expected,b,false,__ATOMIC_RELAXED,__ATOMIC_RELAXED))return;
        }
    }
};
}
