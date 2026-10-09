#include "../src/merge_r8qkeep/family.cpp"
// Frozen B1 oracle from 456fc55: original seen-array and root-shortcut traversal.
#include "task23_family_serial.hpp"
#include <random>
#include <stdexcept>
static void require(bool ok){if(!ok)throw std::runtime_error("family parallel mismatch");}
static void check(const std::vector<std::string>& keys){
    std::vector<int> expected(keys.size());std::iota(expected.begin(),expected.end(),0);
    require(serial_neighbor_families(keys,expected));
    const auto roots=[](std::vector<int> p){std::map<int,int> first;for(std::size_t i=0;i<p.size();++i){int r=i;while(p[r]!=r)r=p[r];first.emplace(r,i);}std::vector<int> out(p.size());for(std::size_t i=0;i<p.size();++i){int r=i;while(p[r]!=r)r=p[r];out[i]=first.at(r);}return out;};
    for(int threads:{1,8})for(std::size_t limit:{0,1,64}){
        std::vector<int> actual(keys.size());std::iota(actual.begin(),actual.end(),0);
        require(neighbor_families(keys,actual,threads,limit));
        require(roots(actual)==roots(expected));
    }
}
int main(){try{
    std::mt19937 rng(2302);const std::string alphabet="ACGT";
    std::vector<std::string> keys;
    for(int family=0;family<32;++family){std::string key(40,'A');for(char& c:key)c=alphabet[rng()%4];
        for(int member=0;member<10;++member){auto k=key;
            for(int m=0;m<member;++m)k[rng()%40]=alphabet[rng()%4];
            if(member%3==0)k[(family*7)%40]='N';if(member%5==0)k[(family*7+13)%40]='N';keys.push_back(k);}
        keys.push_back(key);keys.push_back(key);
    }
    keys.push_back("");check(keys);
    // >64 accepted neighbors, duplicates, N/N mismatches and connected chains.
    std::vector<std::string> dense;std::string seed="ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    for(int p=0;p<40;++p)for(char c:std::string("ACGTN")){auto key=seed;key[p]=c;dense.push_back(key);}
    check(dense);
    std::vector<std::string> chain{seed};for(int group=0;group<5;++group){for(int p=group*8;p<(group+1)*8;++p)seed[p]=seed[p]=='A'?'C':'A';chain.push_back(seed);}check(chain);
    std::vector<int> parent(1);require(!neighbor_families({std::string(40,'X')},parent,8));
    require(!neighbor_families({std::string(40,'N')},parent,8));
    // Check all serialized output and corrected read fields, including edits.
    std::vector<fqdup::r8qkeep::FamilyRead> r1,r2;
    for(std::size_t i=0;i<keys.size();++i){auto key=keys[i];if(key.empty())key="ACGT";
        std::string id="pair"+std::to_string(i),a=key.substr(0,20),b=key.size()>20?key.substr(20):key;
        r1.push_back({"@"+id+"/1",id,a,std::string(a.size(),'5')});r2.push_back({"@"+id+"/2",id,b,std::string(b.size(),'5')});}
    // Three identical high-complexity copies guarantee Q20->Q30 confirmations.
    for(int i=0;i<3;++i){std::string id="confirmed"+std::to_string(i),s="ACGTCAGTGCATGACTCAGT";r1.push_back({"@"+id+"/1",id,s,std::string(20,'5')});r2.push_back({"@"+id+"/2",id,s,std::string(20,'5')});}
    auto a=fqdup::r8qkeep::correct_family(r1,r2,1),b=fqdup::r8qkeep::correct_family(r1,r2,8);
    require(a.edits_tsv==b.edits_tsv&&a.families_tsv==b.families_tsv&&a.stats_json==b.stats_json);
    const auto quiet=fqdup::r8qkeep::correct_family(r1,r2,8,false);
    require(quiet.edits_tsv.empty()&&quiet.families_tsv.empty()&&quiet.stats_json==a.stats_json);
    require(quiet.pairs.size()==a.pairs.size());
    for(std::size_t i=0;i<a.pairs.size();++i)for(int m=0;m<2;++m){const auto& x=a.pairs[i].r[m];const auto& y=quiet.pairs[i].r[m];require(x.h==y.h&&x.id==y.id&&x.s==y.s&&x.q==y.q);}
    require(std::count(a.edits_tsv.begin(),a.edits_tsv.end(),'\n')>1);
    require(a.pairs.size()==b.pairs.size());for(std::size_t i=0;i<a.pairs.size();++i)for(int mate=0;mate<2;++mate){const auto& x=a.pairs[i].r[mate];const auto& y=b.pairs[i].r[mate];require(x.h==y.h&&x.id==y.id&&x.s==y.s&&x.q==y.q);}
    std::cout<<"TASK23_FAMILY_PARALLEL_UNIT_PASS serial-oracle p1/p8 overflow N duplicates chains edits\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
