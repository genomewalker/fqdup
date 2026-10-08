#include "../src/merge_r8qkeep/mixture_fit.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace fqdup::r8qkeep;
static void require(bool condition,const char* message){if(!condition)throw std::runtime_error(message);}
static std::vector<const PairLikelihood*> view(const std::vector<PairLikelihood>& pairs){
    std::vector<const PairLikelihood*> out;for(const auto& p:pairs)out.push_back(&p);return out;
}
static void identical(const MixtureFit& a,const MixtureFit& b){
    require(a.weights==b.weights && a.training_pairs==b.training_pairs &&
        a.unique_likelihoods==b.unique_likelihoods && a.converged==b.converged &&
        a.iterations==b.iterations && a.polish_iterations==b.polish_iterations &&
        a.kkt_gap==b.kkt_gap && a.log_likelihood==b.log_likelihood,"counts order/thread dependence");
}
static void repeated_objective(const std::vector<PairLikelihood>& pairs,const MixtureFit& fit){
    double ll=0;
    std::array<double,mixture_components> gradient{};
    for(const auto& p:pairs){
        const double m=*std::max_element(p.log_likelihood.begin(),p.log_likelihood.end());
        double density=0;
        for(int j=0;j<mixture_components;++j)density+=fit.weights[j]*std::exp(p.log_likelihood[j]-m);
        ll+=m+std::log(density);
        for(int j=0;j<mixture_components;++j)gradient[j]+=std::exp(p.log_likelihood[j]-m)/density;
    }
    require(std::abs(ll-fit.log_likelihood)<1e-9*(1+std::abs(ll)),"weighted objective differs from explicit pairs");
    const double gap=*std::max_element(gradient.begin(),gradient.end())-pairs.size();
    require(std::abs(gap-fit.kkt_gap)<1e-9*pairs.size(),"weighted KKT gradient differs from explicit pairs");
    require(gap<=1.01e-6*pairs.size(),"weighted fit violates KKT tolerance");
}
int main(){
    require(!fit_mixture_counts({}),"empty counts fit");
    std::vector<PairLikelihood> one(19);
    const auto single=fit_mixture_counts(view(one));
    require(single && single->training_pairs==19 && single->unique_likelihoods==1 && single->converged,"single cell count");
    repeated_objective(one,*single);
    const auto single_expanded=fit_mixture_views(view(one));
    for(int j=0;j<mixture_components;++j)require(std::abs(single->weights[j]-single_expanded->weights[j])<1e-14,"single cell expanded weights");

    std::vector<PairLikelihood> pairs;
    for(int kind=0;kind<3;++kind)for(int n=0;n<(kind+1)*11;++n){
        PairLikelihood p;p.id=std::to_string(pairs.size());p.log_likelihood.fill(-20);
        p.log_likelihood[kind]=0;p.log_likelihood[(kind+1)%3]=-2;
        pairs.push_back(std::move(p));
    }
    const auto counts=fit_mixture_counts(view(pairs));
    require(counts && counts->converged && counts->training_pairs==66 && counts->unique_likelihoods==3,"mixed cell counts");
    repeated_objective(pairs,*counts);
    const auto expanded=fit_mixture_views(view(pairs));
    require(expanded && expanded->converged,"expanded fit convergence");
    require(std::abs(expanded->log_likelihood-counts->log_likelihood)<1e-7,"expanded objective mismatch");
    for(int j=0;j<mixture_components;++j)require(std::abs(expanded->weights[j]-counts->weights[j])<2e-5,"expanded weights mismatch");
    require(counts->exp_calls==counts->unique_likelihoods*mixture_components,"counts exp must run once per cell/component");
    // Impossible components, severe scaling, and exact-zero polished weights.
    std::vector<PairLikelihood> sparse(7);
    for(auto& p:sparse){p.log_likelihood.fill(-std::numeric_limits<double>::infinity());p.log_likelihood[0]=-10000;p.log_likelihood[1]=-10800;}
    const auto sparse_fit=fit_mixture_counts(view(sparse));
    require(sparse_fit && sparse_fit->converged && sparse_fit->weights[0]>.999999,"scaled sparse fit");
    repeated_objective(sparse,*sparse_fit);
    identical(*sparse_fit,*fit_mixture_counts(view(sparse),8));
    for(auto& p:sparse)p.log_likelihood.fill(-std::numeric_limits<double>::infinity());
    bool impossible=false;try{fit_mixture_counts(view(sparse));}catch(const std::runtime_error&){impossible=true;}
    require(impossible,"impossible cell accepted");
    identical(*counts,*fit_mixture_counts(view(pairs),8));
    std::reverse(pairs.begin(),pairs.end());
    identical(*counts,*fit_mixture_counts(view(pairs),8));
    for(auto& p:pairs)p.fold=1;
    identical(*counts,*fit_mixture_counts(view(pairs)));
    pairs.front().fold=0;
    bool rejected=false;try{fit_mixture_counts(view(pairs));}catch(const std::invalid_argument&){rejected=true;}
    require(rejected,"mixed folds accepted");

    // Same four base observations, different allocation to whole pairs.
    // A base-count-only statistic would incorrectly identify these datasets.
    LikelihoodAccumulator separated,mixed;
    for(int k=0;k<2;++k){separated.add("a",0,.01,.02,true);separated.add("b",0,.01,.02,false);}
    for(const char* id:{"a","b"}){mixed.add(id,0,.01,.02,true);mixed.add(id,0,.01,.02,false);}
    const auto a=separated.take(),b=mixed.take();
    require(a[0].log_likelihood!=a[1].log_likelihood && b[0].log_likelihood==b[1].log_likelihood,"whole-pair grouping fixture");
    require(fit_mixture_counts(view(a))->unique_likelihoods==2 && fit_mixture_counts(view(b))->unique_likelihoods==1,"base observations were pooled");
    std::cout<<"TASK25_COUNTS_PASS\n";
}
