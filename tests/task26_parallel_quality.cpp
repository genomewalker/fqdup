#include "../src/merge_r8qkeep/compact_store.hpp"
#include <iostream>
using namespace fqdup::r8qkeep;
int main(){
    WorkerPool pool(8);
    for(unsigned alphabet:{1,2,3,5,9}){
        CompactPairStore store;FamilyPair input,output;input.id="widen";
        for(unsigned row=0;row<4096;++row){
            for(unsigned m=0;m<2;++m){auto&r=input.r[m];r.h="@widen/"+std::to_string(m+1);r.id=input.id;r.s.assign(17,'A');r.q.resize(17);
                for(unsigned i=0;i<17;++i)r.q[i]=char(33+(row+m+i)%alphabet);}
            store.append(input.r[0],input.r[1]);
        }
        for(unsigned target:{3,5,9,17}){
            std::string extra;for(unsigned i=0;i<target;++i)extra+=char(33+i);
            store.prepare_qualities(extra,&pool,8);
            CompactPairStore::Decoder decoder;
            for(unsigned row=0;row<4096;++row){store.read(row,output,decoder);for(unsigned m=0;m<2;++m)for(unsigned i=0;i<17;++i)
                if(output.r[m].q[i]!=char(33+(row+m+i)%alphabet))throw std::runtime_error("parallel quality widening");}
        }
    }
    std::cout<<"parallel in-place quality widening all widths PASS\n";
}
