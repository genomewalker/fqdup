#include "../src/merge_r8qkeep/family.cpp"
// Measurement only; cross-check both production key paths.
std::string reason(const Pair& p) {
    if(p.r[0].s.size()<KEYLEN||p.r[1].s.size()<KEYLEN)return "short";
    if(exclude_allg_r2(p.r[1].s))return "allg";
    auto k=p.r[0].s.substr(0,KEYLEN)+p.r[1].s.substr(0,KEYLEN);
    if(std::count(k.begin(),k.end(),'N')>MAXN)return "Ns";
    for(int s=0;s<2;++s){
        const auto& r=p.r[s].s;auto at=r.find(ADAPTERS[s].substr(0,ADSEED));
        if(at!=std::string::npos&&at<KEYLEN)return "adapter";
        auto end=std::min(r.size(),at),g=r.size();
        while(g&&r[g-1]=='G')--g;
        if(r.size()-g>=10)end=std::min(end,g);
        if(low_complexity(std::string_view(r).substr(0,end)))return "complexity";
    }
    return "eligible";
}
int main(int argc,char** argv){
    if(argc!=5)return 2;
    auto a=read_fastq(argv[1]),b=read_fastq(argv[2]);if(a.size()!=b.size())return 3;
    std::vector<Pair> pairs;std::unordered_map<std::string,int> ids;
    for(std::size_t i=0;i<a.size();++i){if(a[i].id!=b[i].id)return 4;ids[a[i].id]=int(i);Pair p;p.r[0]=a[i];p.r[1]=b[i];pairs.push_back(p);}
    auto keys=keys_of<fqdup::r8qkeep::FamilyKey>(pairs);std::vector<std::string> why;
    for(std::size_t i=0;i<pairs.size();++i){auto r=reason(pairs[i]);
        if((r=="eligible")!=(!keys[i].empty())||!(keys[i]==fqdup::r8qkeep::compact_family_key(pairs[i])))return 5;
        why.push_back(r);
    }
    std::unordered_map<std::string,std::vector<int>> truth;
    std::ifstream tf(argv[3]);std::string id,f,line;
    while(tf>>id>>f){if(!ids.count(id))return 6;truth[f].push_back(ids.at(id));}
    std::vector<int> group(pairs.size(),-1);std::ifstream gf(argv[4]);std::getline(gf,line);int gi=0;
    while(std::getline(gf,line)){auto at=line.find_last_of('\t');std::istringstream members(line.substr(at+1));while(std::getline(members,id,',')){if(!ids.count(id))return 7;group[ids.at(id)]=gi;}++gi;}
    std::cout<<"id1\tid2\treason1\treason2\n";
    for(const auto& entry:truth){const auto& m=entry.second;if(m.size()<2||m.size()>50)continue;
        for(std::size_t i=0;i<m.size();++i)for(std::size_t j=i+1;j<m.size();++j){int x=m[i],y=m[j];
            if(group[x]>=0&&group[x]==group[y])continue;
            auto kx=pairs[x].r[0].s.substr(0,20)+pairs[x].r[1].s.substr(0,20);
            auto ky=pairs[y].r[0].s.substr(0,20)+pairs[y].r[1].s.substr(0,20);
            if(kx==ky&&kx.find('N')==std::string::npos)std::cout<<pairs[x].r[0].id<<'\t'<<pairs[y].r[0].id<<'\t'<<why[x]<<'\t'<<why[y]<<'\n';
        }
    }
}
