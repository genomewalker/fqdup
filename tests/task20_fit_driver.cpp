#include "../src/merge_r8qkeep/mixture_fit.hpp"
#include <fstream>
#include <iostream>
#include <stdexcept>
int main(int argc,char** argv) {
    try {
        if(argc!=3&&argc!=4)throw std::runtime_error("usage: fit_driver CENSUS.tsv NEW_PRIOR.txt [THREADS]");
        const int threads=argc==4?std::stoi(argv[3]):1;
        std::ifstream input(argv[1]);if(!input)throw std::runtime_error("cannot open census");
        auto pairs=fqdup::r8qkeep::read_census_likelihood(input);
        std::array<std::optional<fqdup::r8qkeep::MixtureFit>,2> fits;
        for(int f=0;f<2;++f){std::vector<fqdup::r8qkeep::PairLikelihood> fold;for(const auto& p:pairs)if(p.fold==f)fold.push_back(p);fits[f]=fqdup::r8qkeep::fit_mixture(fold,threads);if(fits[f])std::cerr<<"FIT_CACHE fold="<<f<<" pairs="<<fits[f]->training_pairs<<" keys="<<fits[f]->unique_likelihoods<<" calls="<<fits[f]->exp_calls<<" hits="<<fits[f]->exp_cache_hits<<'\n';}
        std::ofstream output(argv[2]);output<<fqdup::r8qkeep::prior_text(fits);if(!output)throw std::runtime_error("prior write failed");
        return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
