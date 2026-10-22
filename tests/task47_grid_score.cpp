// Same truth conventions as task47_attribution, all cutoffs in one truth pass.
#define main task47_original_attribution_main
#include "task47_attribution.cpp"
#undef main
#include <algorithm>
#include <set>
char complement(char b){return b=='A'?'T':b=='T'?'A':b=='C'?'G':b=='G'?'C':b;}
struct GridMetrics {
    MateStats mates;
    Channel channel;
    long long merged=0,fp=0,wrong_len=0,exact=0,base_err=0,fn=0;
    long long first10=0,last10=0,c5=0,g3=0;
};
int main(int argc,char** argv){try{
    require(argc==5,"usage: grid_score truth.jsonl output-root cutoffs.txt library");
    const string root=argv[2],lib=argv[4];ifstream list(argv[3]);require(bool(list),"grid list");
    vector<Arm> arms;string cut;while(getline(list,cut)){if(cut.empty())continue;arms.emplace_back();arms.back().name=cut;arms.back().tagged=true;load(arms.back(),root+"/"+cut,"","");}
    require(!arms.empty(),"empty grid");
    unordered_map<string,array<map<size_t,Edit>,3>> family;
    ifstream ledger(root+"/debug/family-edits.tsv");require(bool(ledger),"family ledger");string line;getline(ledger,line);
    while(getline(ledger,line)){istringstream row(line);string id;int mate,oq,q;size_t p;char before,after;require(bool(row>>id>>mate>>p>>before>>oq>>after>>q),"family row");if(before!=after)family[id][mate].emplace(p,Edit{before,after});}
    vector<GridMetrics> metrics(arms.size());ifstream truth(argv[1]);require(bool(truth),"truth open");set<string> seen;
    while(getline(truth,line)){
        Parser parser{line};auto j=parser.value();const string id=j.at("id").s,target=j.at("merged_target").s;
        require(seen.insert(id).second,"duplicate truth id");const bool should=target!="null"&&target!="-"&&!target.empty();
        const long long length=stoll(j.at("frag_len").s);const auto& reads=j.at("reads").a;
        map<long long,char> opportunities;
        for(int mate=0;mate<2;++mate)for(const auto& o:reads[mate].at("origins").a)if(o.at("kind").s=="molecule")opportunities[stoll(o.at("coordinate").s)]=mate?complement(o.at("truth_base").s[0]):o.at("truth_base").s[0];
        unordered_map<long long,size_t> index;vector<long long> coordinates;vector<char> bases;
        for(const auto& o:opportunities){index[o.first]=coordinates.size();coordinates.push_back(o.first);bases.push_back(o.second);}
        vector<array<long long,2>> copies(coordinates.size());
        for(size_t ai=0;ai<arms.size();++ai){
            auto& arm=arms[ai];auto& m=metrics[ai];fill(copies.begin(),copies.end(),array<long long,2>{0,0});
            const auto addcopy=[&](long long coord,bool good){auto it=index.find(coord);require(it!=index.end(),"copy outside opportunity");++copies[it->second][good?0:1];};
            auto merged=arm.merged.find(id);auto mates=arm.mates.find(id);
            if(merged!=arm.merged.end()){
                require(mates==arm.mates.end(),"merged/mate collision");++m.merged;
                const auto& sequence=merged->second;
                if(!should){++m.fp;m.channel.synthetic+=sequence.size();}
                else if(sequence.size()!=target.size()){++m.wrong_len;m.channel.synthetic+=sequence.size();}
                else{bool wrong=false;for(size_t p=0;p<sequence.size();++p){const bool good=sequence[p]==target[p];wrong|=!good;addcopy(p,good);}wrong?++m.base_err:++m.exact;}
                arm.merged.erase(merged);
            }else{
                m.fn+=should;
                if(mates!=arm.mates.end()){
                    for(const auto& r:mates->second){
                        require(r.mate==1||r.mate==2,"unidentified mate");const auto& raw=reads[r.mate-1];string sequence=raw.at("sequence").s;const auto& origins=raw.at("origins").a;
                        auto f=family.find(id);if(f!=family.end())for(const auto& e:f->second[r.mate]){require(sequence.at(e.first)==e.second.before,"family old base");sequence[e.first]=e.second.after;}
                        const auto tag=r.h.rfind(" FQT=v1;");require(tag!=string::npos,"missing FQT");const auto span=fields(r.h.substr(tag+1,r.h.find(' ',tag+1)-tag-1));
                        const auto begin=stoull(span.at("s")),end=stoull(span.at("e"));require(end>=begin&&end-begin==r.s.size()&&end<=sequence.size(),"span invalid");
                        map<size_t,char> edits;const auto qc=r.h.rfind(" FQC=v1;edits=");
                        if(qc!=string::npos){istringstream input(r.h.substr(qc+14));string edit;while(getline(input,edit,';')){replace(edit.begin(),edit.end(),':',' ');istringstream row(edit);size_t p,donor;char base;int q;require(bool(row>>p>>donor>>base>>q),"FQC syntax");require(p>=begin&&p<end&&r.q[p-begin]==char(q+33),"FQC range/quality");edits.emplace(p,base);}}
                        ++m.mates.n;m.mates.raw+=sequence.size();bool bad=false;
                        for(size_t p=0;p<sequence.size();++p){const auto& o=origins[p];const bool molecule=o.at("kind").s=="molecule";require(molecule||o.at("kind").s=="adapter","unknown origin");
                            if(p>=begin&&p<end){const char call=r.s[p-begin];auto edit=edits.find(p);require(call==(edit==edits.end()?sequence[p]:edit->second),"untraced base change");
                                if(molecule){const bool good=call==o.at("truth_base").s[0];good?++m.mates.k:++m.mates.em;bad|=!good;addcopy(stoll(o.at("coordinate").s),good);}else{++m.mates.ea;bad=true;}}
                            else if(molecule)++m.mates.t;else ++m.mates.at;
                        }
                        m.mates.bad+=bad;
                    }
                    arm.mates.erase(mates);
                }
            }
            m.channel.opportunity+=coordinates.size();m.channel.adapter=m.mates.ea;
            for(size_t p=0;p<copies.size();++p){const auto good=copies[p][0],bad=copies[p][1];m.channel.copies_ok+=good;m.channel.copies_bad+=bad;
                if(good&&bad)++m.channel.cov_conflict;else if(good)++m.channel.cov_ok;else if(bad)++m.channel.cov_bad;else ++m.channel.omitted;
                if(good||bad){const auto c=coordinates[p];m.first10+=c>=0&&c<10;m.last10+=c>=std::max(0LL,length-10)&&c<length;m.c5+=c>=0&&c<10&&bases[p]=='C';m.g3+=c>=std::max(0LL,length-10)&&c<length&&bases[p]=='G';}
            }
        }
    }
    cout<<"library\tcutoff\tmates\tK\tEm\tEa\tT\tRbad\tadapter_trimmed\traw\tmerged\tFP\twrong_len\texact\tbase_err\tFN\tsynthetic\tcovered_correct\tcovered_wrong\tconflicting\tfirst10\tlast10\tC5\tG3\n";
    for(size_t i=0;i<arms.size();++i){require(arms[i].merged.empty()&&arms[i].mates.empty(),"output IDs absent from truth");const auto& m=metrics[i];require(m.mates.k+m.mates.em+m.mates.ea+m.mates.t+m.mates.at==m.mates.raw,"mate conservation");
        cout<<lib<<'\t'<<arms[i].name<<'\t'<<m.mates.n<<'\t'<<m.mates.k<<'\t'<<m.mates.em<<'\t'<<m.mates.ea<<'\t'<<m.mates.t<<'\t'<<m.mates.bad<<'\t'<<m.mates.at<<'\t'<<m.mates.raw<<'\t'<<m.merged<<'\t'<<m.fp<<'\t'<<m.wrong_len<<'\t'<<m.exact<<'\t'<<m.base_err<<'\t'<<m.fn<<'\t'<<m.channel.synthetic<<'\t'<<m.channel.cov_ok<<'\t'<<m.channel.cov_bad<<'\t'<<m.channel.cov_conflict<<'\t'<<m.first10<<'\t'<<m.last10<<'\t'<<m.c5<<'\t'<<m.g3<<'\n';}
}catch(const exception& e){cerr<<e.what()<<'\n';return 1;}}
