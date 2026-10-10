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
    // Bulk corrections must preserve logical-to-physical mapping and the last
    // value at duplicate positions, including edits to existing exceptions.
    struct Edit {unsigned pair,position,mate;char base,quality;bool base_changed=true;};
    for(bool canonical:{false,true})for(int workers:{1,8}){
        CompactPairStore bulk,scalar;
        for(unsigned row=0;row<4096;++row){
            const auto id=std::to_string(4096-row);std::string seq(137,'A'),q(137,'I');seq[7]='N';
            bulk.append("@"+id+"/1 comment",seq,q,"@"+id+"/2 comment",seq,q);
            scalar.append("@"+id+"/1 comment",seq,q,"@"+id+"/2 comment",seq,q);
        }
        if(canonical){bulk.canonicalize();scalar.canonicalize();}
        auto generic=bulk;
        std::vector<Edit> edits;
        for(unsigned i=0;i<1200;++i){const auto row=(i*239)%4096;edits.push_back({row,7,i%2,'T',char(33+i%70)});edits.push_back({row,7,i%2,'T',char(33+(i+3)%70),false});edits.push_back({row,7,i%2,'C',char(33+(i+1)%70)});edits.push_back({row,136,(i+1)%2,'N','!'});}
        for(const auto& e:edits)scalar.edit(e.pair,e.mate,e.position,e.base,e.quality);
        generic.apply_edits(edits,&pool,workers);
        bulk.apply_edits(std::move(edits),&pool,workers,[](const Edit& e){return e.base_changed;});
        FamilyPair a,b,c;CompactPairStore::Decoder da,db,dc;
        for(unsigned row=0;row<4096;++row){bulk.read(row,a,da);scalar.read(row,b,db);generic.read(row,c,dc);for(unsigned m=0;m<2;++m)if(c.r[m].s!=a.r[m].s||c.r[m].q!=a.r[m].q)throw std::runtime_error("filtered versus generic edits");if(a.id!=b.id)throw std::runtime_error("bulk ID");
            for(unsigned m=0;m<2;++m)if(a.r[m].h!=b.r[m].h||a.r[m].s!=b.r[m].s||a.r[m].q!=b.r[m].q)throw std::runtime_error("bulk corrected content");
            bulk.content(row,b,db);for(unsigned m=0;m<2;++m)if(a.r[m].s!=b.r[m].s||a.r[m].q!=b.r[m].q)throw std::runtime_error("content-only decoder");}
    }
    std::cout<<"parallel block quality edits and alphabet changes PASS\n";
}
