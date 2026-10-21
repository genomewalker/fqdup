#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <set>
#include <vector>
struct Family {
    std::string id, size, disbanded;
    unsigned merged[2]{};
    std::set<unsigned> lengths[2];
};
int main(int argc,char** argv) {
    try {
        if(argc!=4)throw std::runtime_error("parent merged.fq, candidate merged.fq, families.tsv required");
        std::ifstream input(argv[3]);
        if(!input)throw std::runtime_error("open families");
        std::vector<Family> families;
        std::unordered_map<std::string,std::size_t> members;
        std::string line;
        std::getline(input,line);
        while(std::getline(input,line)) {
            std::istringstream row(line);Family f;std::string diameter,ids;
            if(!(row>>f.id>>f.size>>diameter>>f.disbanded>>ids))throw std::runtime_error("invalid family row");
            std::istringstream list(ids);std::string id;
            while(std::getline(list,id,','))members.emplace(id,families.size());
            families.push_back(std::move(f));
        }
        for(int arm=0;arm<2;++arm) {
            std::ifstream fq(argv[arm+1]);
            if(!fq)throw std::runtime_error("open merged FASTQ");
            std::string h,s,p,q;
            while(std::getline(fq,h)) {
                if(!std::getline(fq,s)||!std::getline(fq,p)||!std::getline(fq,q))throw std::runtime_error("truncated FASTQ");
                auto id=h.substr(1,h.find_first_of(" \t")-1);
                if(id.size()>2&&id[id.size()-2]=='/')id.resize(id.size()-2);
                auto it=members.find(id);if(it==members.end())continue;
                auto& f=families[it->second];++f.merged[arm];f.lengths[arm].insert(s.size());
            }
        }
        for(const auto& f:families)for(int arm=0;arm<2;++arm)
            std::cout<<f.id<<'\t'<<f.size<<'\t'<<f.disbanded<<'\t'<<arm+1<<'\t'<<f.merged[arm]<<'\t'<<f.lengths[arm].size()<<'\t'
                <<(f.merged[arm]>0&&f.merged[arm]<std::stoul(f.size)?"split":"concordant")<<'\n';
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
