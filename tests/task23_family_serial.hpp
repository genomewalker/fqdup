bool serial_neighbor_families(const std::vector<std::string>& keys,std::vector<int>& parent){
    const auto code=[](char c){switch(c){case 'A':return 0;case 'C':return 1;case 'G':return 2;case 'T':return 3;case 'N':return 4;default:return -1;}};
    for(const auto& key:keys)if(!key.empty()){if(key.size()!=40)return false;for(char c:key)if(code(c)<0)return false;}
    const auto started=std::chrono::steady_clock::now();
    const auto root=[&](int i){while(parent[i]!=i){parent[i]=parent[parent[i]];i=parent[i];}return i;};
    const auto unite=[&](int i,int j){i=root(i);j=root(j);if(i!=j)parent[i]=j;};
    std::unordered_map<std::string,int> exact;std::vector<int> representative;
    for(std::size_t i=0;i<keys.size();++i)if(!keys[i].empty()){auto [it,new_key]=exact.emplace(keys[i],representative.size());if(new_key)representative.push_back(i);else unite(i,representative[it->second]);}
    using Index=std::unordered_map<std::uint64_t,std::vector<int>>;
    std::array<Index,3> index;const int bounds[]={0,13,26,40};
    const auto packed=[&](const std::string& key,int part){std::uint64_t x=0;for(int p=bounds[part];p<bounds[part+1];++p)x|=std::uint64_t(code(key[p]))<<(3*(p-bounds[part]));return x;};
    std::size_t cells=0,max_bucket=0;
    for(int part=0;part<3;++part){for(std::size_t i=0;i<representative.size();++i)index[part][packed(keys[representative[i]],part)].push_back(i);cells+=index[part].size();for(const auto& entry:index[part])max_bucket=std::max(max_bucket,entry.second.size());}
    const auto generated=std::chrono::steady_clock::now();
    std::vector<std::size_t> seen(representative.size(),0);std::uint64_t lookups=0,candidates=0,verified=0;
    for(std::size_t query=0;query<representative.size();++query){const auto& key=keys[representative[query]];
        for(int part=0;part<3;++part){const int length=bounds[part+1]-bounds[part];std::vector<int> normal,npos;
            for(int p=0;p<length;++p)(key[bounds[part]+p]=='N'?npos:normal).push_back(p);
            std::vector<std::uint64_t> bases{packed(key,part)};
            for(int p:npos){std::vector<std::uint64_t> next;for(auto base:bases)for(int c=0;c<5;++c)next.push_back((base&~(std::uint64_t(7)<<(3*p)))|(std::uint64_t(c)<<(3*p)));bases=std::move(next);}
            const int budget=2-int(npos.size());
            const auto lookup=[&](std::uint64_t value){++lookups;const auto it=index[part].find(value);if(it==index[part].end())return;candidates+=it->second.size();
                for(int other:it->second){if(std::size_t(other)>=query||seen[other]==query+1)continue;seen[other]=query+1;++verified;
                    if(root(representative[query])==root(representative[other]))continue;
                    if(mismatches(key,keys[representative[other]])<=MAXMIS)unite(representative[query],representative[other]);}};
            for(auto base:bases){lookup(base);if(budget<1)continue;
                for(std::size_t a=0;a<normal.size();++a){const int p=normal[a],original=code(key[bounds[part]+p]);
                    for(int c=0;c<5;++c)if(c!=original){const auto first=(base&~(std::uint64_t(7)<<(3*p)))|(std::uint64_t(c)<<(3*p));lookup(first);if(budget<2)continue;
                        for(std::size_t b=a+1;b<normal.size();++b){const int q=normal[b],second=code(key[bounds[part]+q]);for(int d=0;d<5;++d)if(d!=second)lookup((first&~(std::uint64_t(7)<<(3*q)))|(std::uint64_t(d)<<(3*q)));}}}
            }
        }
    }
    if(std::getenv("FQDUP_FAMILY_DIAGNOSTICS")){const auto finished=std::chrono::steady_clock::now();std::cerr<<"FAMILY_INDEX pairs="<<keys.size()<<" unique_keys="<<representative.size()<<" cells="<<cells<<" max_bucket="<<max_bucket<<" lookups="<<lookups<<" candidates="<<candidates<<" verified="<<verified<<" generation_seconds="<<std::chrono::duration<double>(generated-started).count()<<" verification_seconds="<<std::chrono::duration<double>(finished-generated).count()<<'\n';}
    return true;
}

