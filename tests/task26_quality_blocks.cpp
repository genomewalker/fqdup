#include "../src/merge_r8qkeep/quality_blocks.hpp"
#include <iostream>
#include <random>
using namespace fqdup::r8qkeep;
static void require(bool ok){if(!ok)throw std::runtime_error("quality block mismatch");}
static void check(QualityBlocks& blocks,const std::string& expected,QualityBlocks::Cursor& cursor){
    require(blocks.size()==expected.size());std::string out;
    for(std::size_t at=0;at<expected.size();at+=151){const auto n=std::min<std::size_t>(151,expected.size()-at);blocks.read(at,n,out,cursor);require(out==expected.substr(at,n));}
}
int main(){
    std::mt19937 rng(2605);WorkerPool pool(8);
    for(unsigned alphabet:{1,2,3,5,9,16,17,94,256})for(bool runs:{false,true}){
        std::string expected;
        for(unsigned i=0;i<(alphabet==17&&!runs?257:9)*QualityBlocks::block_size+513;++i)expected+=char((runs?i/137:rng())%alphabet);
        QualityBlocks blocks;std::size_t at=0;
        while(at<expected.size()){const auto n=std::min<std::size_t>(expected.size()-at,1+rng()%512);blocks.append(std::string_view(expected).substr(at,n));at+=n;}
        QualityBlocks::Cursor cursor;check(blocks,expected,cursor);blocks.finish();check(blocks,expected,cursor);
        std::string out;
        for(unsigned length=0;length<=512;++length){const auto first=rng()%(expected.size()-length+1);blocks.read(first,length,out,cursor);require(out==expected.substr(first,length));}
        for(unsigned i=0;i<300;++i){const auto pos=rng()%expected.size();const auto value=char(rng()%256);blocks.edit(pos,static_cast<unsigned char>(value));expected[pos]=value;
            const auto first=pos>128?pos-128:0;blocks.read(first,std::min<std::size_t>(512,expected.size()-first),out,cursor);require(out==expected.substr(first,out.size()));}
        check(blocks,expected,cursor);
        for(int workers:{1,8}){
            std::vector<QualityBlocks::Edit> edits;
            for(unsigned i=0;i<3000;++i){const auto pos=rng()%expected.size();const auto value=static_cast<unsigned char>(rng());edits.push_back({pos,value});expected[pos]=char(value);
                if(i%17==0){edits.push_back({pos,static_cast<unsigned char>(value^128)});expected[pos]=char(value^128);}}
            blocks.rewrite(std::move(edits),workers==1?nullptr:&pool,workers);check(blocks,expected,cursor);
        }
        QualityBlocks joined;std::string prefix="prefix";joined.append(prefix);joined.join(blocks);prefix+=expected;check(joined,prefix,cursor);
        QualityBlocks extra;extra.append("pending only");joined.join(extra);prefix+="pending only";check(joined,prefix,cursor);
    }
    QualityBlocks small;std::string text;
    for(unsigned i=0;i<500;++i){QualityBlocks part;const std::string fragment(i%37+1,char(i%256));part.append(fragment);if(i%3==0)part.finish();small.join(part);text+=fragment;}
    QualityBlocks::Cursor cursor;check(small,text,cursor);
    QualityBlocks constant;constant.append(std::string(100000,'I'));constant.finish();require(constant.owned_bytes()<100000);
    std::cout<<"quality blocks: alphabets, all lengths0..512, runs/raw, joins, edits/cache invalidation and p1/p8 rewrite PASS\n";
}
