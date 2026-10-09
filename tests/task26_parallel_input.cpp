#include "../src/merge_r8qkeep/compact_store.hpp"
#include "fqdup/fastq_types.hpp"
#include <atomic>
#include <iostream>
using namespace fqdup::r8qkeep;
class MockReader:public FastqReaderBase {
    std::size_t at_=0,count_;unsigned mate_;bool fail_;std::atomic<bool> reading_{false};
public:
    MockReader(std::size_t count,unsigned mate,bool fail):count_(count),mate_(mate),fail_(fail){}
    bool read(FastqRecord& r)override{
        if(reading_.exchange(true))throw std::runtime_error("reader used concurrently");
        struct Reset{std::atomic<bool>& flag;~Reset(){flag=false;}}reset{reading_};
        if(fail_&&at_==4100)throw std::runtime_error("injected reader failure");
        if(at_==count_)return false;
        r.header="@raw-"+std::to_string(at_)+"/"+std::to_string(mate_)+" comment";r.seq.clear();r.qual.clear();r.plus="+";
        for(std::size_t i=0;i<at_%513;++i){r.seq+="ACGTN"[(at_+i+mate_)%5];r.qual+=char(33+(at_+i+mate_)%22);}++at_;return true;
    }
    std::uint64_t record_count()const override{return at_;}
    const char* backend_name()const override{return "test";}
};
std::unique_ptr<FastqReaderBase> make_chained_fastq_reader(const std::vector<std::string>& paths,std::size_t){
    return std::make_unique<MockReader>(std::stoul(paths.at(0)),std::stoul(paths.at(1)),paths.size()>2);
}
namespace {
#include "../src/merge_r8qkeep/parallel_input.inc"
}
int main(){
    for(int threads:{2,8}){
        auto store=parallel_compact_input({"12303","1"},{"12303","2"},threads);
        if(store->size()!=12303)throw std::runtime_error("parallel cardinality");
        MockReader a(12303,1,false),b(12303,2,false);FastqRecord x,y;FamilyPair observed;CompactPairStore::Decoder cursor;
        for(std::size_t i=0;i<store->size();++i){a.read(x);b.read(y);store->read(i,observed,cursor);
            if(observed.r[0].h!=x.header||observed.r[0].s!=x.seq||observed.r[0].q!=x.qual||observed.r[1].h!=y.header||observed.r[1].s!=y.seq||observed.r[1].q!=y.qual)throw std::runtime_error("parallel input parity");}
    }
    for(bool fail:{false,true}){bool caught=false;try{
        if(fail)parallel_compact_input({"12303","1"},{"12303","2","fail"},8);
        else parallel_compact_input({"12303","1"},{"12304","2"},8);
    }catch(const std::runtime_error& e){caught=std::string(e.what())==(fail?"injected reader failure":"pair count");}if(!caught)throw std::runtime_error("reader failure missing");}
    auto empty=parallel_compact_input({"0","1"},{"0","2"},8);if(!empty->empty())throw std::runtime_error("empty input");
    std::cout<<"parallel input parity, cardinality, empty input and failure drain PASS\n";
}
