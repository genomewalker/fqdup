#include "../src/merge_r8qkeep/family.cpp"
// Frozen B1 oracle from 456fc55: original seen-array and root-shortcut traversal.
#include "task23_family_serial.hpp"
#include "task25_columns_ref.hpp"
#include "task25_family_counters_ref.hpp"
#include "task25_family_keys_ref.hpp"
#include <cfenv>
#include <cstring>
#include <random>
#include <stdexcept>
static void require(bool ok){if(!ok)throw std::runtime_error("family parallel mismatch");}
static void packed_key_parity(){
    std::mt19937 rng(2560);const int bounds[]{0,13,26,40};
    for(int trial=0;trial<10000;++trial){
        std::string a(40,'A'),b(40,'A');
        for(char& c:a)c="ACGTN"[rng()%5];for(char& c:b)c="ACGTN"[rng()%5];
        const PackedFamilyKey x(a),y(b);
        for(int p=0;p<3;++p)require(x.part(p)==fqdup::r8qkeep::FamilyPartIndex::pack(a,bounds[p],bounds[p+1]));
        const auto lo=PackedFamilyKey::differences(x.lo,y.lo),hi=PackedFamilyKey::differences(x.hi,y.hi);
        require(__builtin_popcountll(lo)+__builtin_popcountll(hi)==mismatches(a,b));
        require(__builtin_popcountll(lo&((1ULL<<39)-1))==mismatches(a.substr(0,13),b.substr(0,13)));
        require(__builtin_popcountll(lo>>39)+__builtin_popcountll(hi&((1ULL<<18)-1))==mismatches(a.substr(13,13),b.substr(13,13)));
    }
}
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
static void column_parity(){
    std::mt19937 rng(2516);
    for(int size=0;size<=MAXCOMP;++size)for(int trial=0;trial<32;++trial){
        std::vector<Pair> pairs(24);
        for(std::size_t id=0;id<pairs.size();++id)for(int side=0;side<2;++side){
            auto& read=pairs[id].r[side];read.s.resize(1+rng()%30);
            for(auto& base:read.s)base="ACGTN"[rng()%5];
            read.q.resize(read.s.size());for(char& q:read.q)q=char(33+rng()%61);
        }
        std::vector<int> group(24);std::iota(group.begin(),group.end(),0);
        std::shuffle(group.begin(),group.end(),rng);group.resize(size);
        // Exercise all bounded base slots, not only the usual four DNA bases.
        if(trial==0)for(int member=0;member<size;++member)
            for(int side=0;side<2;++side)pairs[group[member]].r[side].s[0]=char('A'+member);
        std::array<const std::vector<Column>*,2> prepared{};
        require(second_molecule(pairs,group,&prepared)==frozen_family_counters::second_molecule(pairs,group));
        require(high_quality_counts(pairs,group,prepared)==frozen_family_counters::high_quality_counts(pairs,group));
        for(int side=0;side<2;++side){
            const auto expected=frozen_columns::columns(pairs,group,side);
            const auto& actual=size>=2?*prepared[side]:columns(pairs,group,side);
            require(expected.size()==actual.size());
            for(std::size_t c=0;c<actual.size();++c){
                require(expected[c].n==actual[c].n&&expected[c].top==actual[c].top&&
                        expected[c].votes.size()==actual[c].votes.size());
                for(std::size_t v=0;v<actual[c].votes.size();++v){
                    const auto vote=actual[c].votes[v];
                    require(vote.first==expected[c].votes[v].first);
                    std::vector<int> ids;for(int id:vote.second)ids.push_back(id);
                    require(ids==expected[c].votes[v].second&&ids.size()==vote.second.size());
                    for(std::size_t i=0;i<ids.size();++i)require(ids[i]==vote.second[i]);
                }
            }
        }
    }
    // Keep disagreement below MAXDIS so shared minority pairs, not the
    // preliminary rate check, decide whether these shuffled-ID groups disband.
    for(int pattern=0;pattern<3;++pattern){
        std::vector<Pair> pairs(24);const std::vector<int> group{23,7,18,2,9,0};
        for(auto& pair:pairs)for(auto& read:pair.r){read.s.assign(50,'A');read.q.assign(50,'I');}
        for(int member:{1,4})pairs[group[member]].r[0].s[3]='C';
        if(pattern)for(int member:{1,pattern==1?4:5})pairs[group[member]].r[1].s[8]='C';
        require(second_molecule(pairs,group)==(pattern==1));
        require(second_molecule(pairs,group)==frozen_family_counters::second_molecule(pairs,group));
        std::array<const std::vector<Column>*,2> prepared{};
        require(high_quality_counts(pairs,group,prepared)==frozen_family_counters::high_quality_counts(pairs,group));
    }
    std::cout<<"FAMILY_COLUMN_BYTES compact="<<sizeof(Column)
             <<" old_object="<<sizeof(frozen_columns::Column)
             <<" old_nested_allocations=excluded"<<std::endl;
}
static void key_parity(){
    std::mt19937 rng(2538);std::vector<Pair> all(1025);
    const std::array<std::size_t,10> lengths{0,1,19,20,21,31,64,100,150,512};
    for(std::size_t id=0;id<all.size();++id)for(int side=0;side<2;++side){
        auto& sequence=all[id].r[side].s;sequence.resize(lengths[(id+side*3)%lengths.size()]);
        for(char& c:sequence)c="ACGT"[rng()%4];
        switch((id/lengths.size())%8){
        case 1:for(std::size_t k=0;k<std::min<std::size_t>(id%4,sequence.size());++k)sequence[k]='N';break;
        case 2:sequence.assign(sequence.size(),'A');break;
        case 3:{const std::size_t pos=std::array<std::size_t,5>{0,19,20,29,40}[id%5];if(pos+ADSEED<=sequence.size())sequence.replace(pos,ADSEED,ADAPTERS[side].substr(0,ADSEED));break;}
        case 4:{const auto tail=std::min<std::size_t>(9+id%3,sequence.size());sequence.replace(sequence.size()-tail,tail,tail,'G');break;}
        case 5:for(std::size_t k=0;k<sequence.size();++k)sequence[k]="ACGT"[k%4];break;
        case 6:if(!sequence.empty())sequence[id%sequence.size()]='X';break;
        default:break;
        }
    }
    const int rounding=std::fegetround();std::size_t cases=0;
    for(int mode:{FE_TONEAREST,FE_DOWNWARD,FE_UPWARD,FE_TOWARDZERO}){
        require(std::fesetround(mode)==0);
        for(int workers:{2,8}){
            fqdup::r8qkeep::WorkerPool pool(workers);
            for(std::size_t count:{0,1,127,128,129,511,512,513,1023,1024,1025}){
                std::vector<Pair> pairs(all.begin(),all.begin()+count);
                const auto expected=frozen_keys_of(pairs);
                require(keys_of(pairs)==expected);
                require(keys_of(pairs,&pool,workers)==expected);++cases;
            }
        }
    }
    require(std::fesetround(rounding)==0);
    std::cout<<"TASK25_FAMILY_KEYS_PASS cases="<<cases<<'\n';
}

int main(){
    packed_key_parity();try{
    // Independent complete-window oracle for allocation-free rolling DUST.
    std::mt19937 dust_rng(2547);
    for(int trial=0;trial<200;++trial){
        std::string seq(trial*7%513,'A');
        if(trial%4)for(char& c:seq)c="ACGTN?"[dust_rng()%(trial%4==1?4:6)];
        bool expected=false;
        if(seq.size()>=5){const int nt=int(seq.size())-2,wt=std::min(62,nt);
            for(int start=0;start+wt<=nt;++start){std::array<int,64> counts{};long valid=0,score=0;
                for(int i=start;i<start+wt;++i){int word=0;bool good=true;
                    for(int k=0;k<3;++k){const auto p=std::string_view("ACGT").find(seq[i+k]);if(p==std::string_view::npos){good=false;break;}word=word*4+int(p);}
                    if(good){++counts[word];++valid;}
                }
                for(int n:counts)score+=n*(n-1)/2;
                expected|=valid>=2&&double(score)/(double(valid)*(valid-1)/2)>DUST_MAX;
            }
        }
        require(low_complexity(seq)==expected);
    }
    key_parity();
    const int rounding=std::fegetround();
    for(int mode:{FE_TONEAREST,FE_DOWNWARD,FE_UPWARD,FE_TOWARDZERO}){
        require(std::fesetround(mode)==0);const auto errors=family_quality_errors();
        Read read;read.q.resize(1);
        for(int byte=0;byte<256;++byte){read.q[0]=char(byte);double expected=std::pow(10.0,-qv(read,0)/10.0);
            require(std::memcmp(&expected,&errors[byte],sizeof(double))==0);}
    }
    require(std::fesetround(rounding)==0);
    column_parity();
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
