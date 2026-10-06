#include "mixture_fit.hpp"
#include "shortest_double.hpp"
#include "fdlibm.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
namespace fqdup::r8qkeep {
namespace {
using Values = std::array<double,mixture_components>;
double sum(const Values& values) { double out=0;for(double v:values)out+=v;return out; }
double logsum(const Values& values) {
    const double m=*std::max_element(values.begin(),values.end());
    if(!std::isfinite(m))throw std::runtime_error("all mixture components impossible");
    double s=0;for(double x:values)s+=fqdup::fdlibm::exp(x-m);
    return m+fqdup::fdlibm::log(s);
}
}
std::vector<PairLikelihood> read_census_likelihood(std::istream& census) {
    const auto likelihood=[](double d,bool agree){return fqdup::fdlibm::log(std::max(1e-15,std::min(1-1e-15,agree?1-d:d)));};
    std::array<std::array<double,2>,mixture_components> flat{};
    for(int j=1;j<mixture_components;++j){const double r=double(j-1)/64,e=.75*(r*r),d=2*e-4*e*e/3;flat[j]={likelihood(d,true),likelihood(d,false)};}
    std::vector<PairLikelihood> pairs;std::unordered_map<std::string,std::size_t> slots;
    std::string line;if(!std::getline(census,line))throw std::runtime_error("missing census header");
    while(std::getline(census,line)) {
        if(line.empty())continue;
        std::istringstream row(line);std::vector<std::string> fields;std::string field;while(std::getline(row,field,'\t'))fields.push_back(field);
        if(fields.size()!=11)throw std::runtime_error("census column count");
        const int fold=std::stoi(fields[1]);if(fold!=0&&fold!=1)throw std::runtime_error("census fold");
        auto [it,inserted]=slots.emplace(fields[0],pairs.size());if(inserted){PairLikelihood p;p.id=fields[0];p.fold=fold;pairs.push_back(std::move(p));}
        auto& p=pairs[it->second];if(p.fold!=fold)throw std::runtime_error("fold mismatch");
        const double e1=std::stod(fields[6]),e2=std::stod(fields[7]);const bool agree=fields[8]==fields[9];
        p.log_likelihood[0]+=likelihood(e1+e2-4*e1*e2/3,agree);
        for(int j=1;j<mixture_components;++j)p.log_likelihood[j]+=flat[j][agree?0:1];
    }
    if(census.bad())throw std::runtime_error("census read failure");
    return pairs;
}
std::optional<MixtureFit> fit_mixture(const std::vector<PairLikelihood>& pairs) {
    if(pairs.empty())return std::nullopt;
    MixtureFit result;result.training_pairs=pairs.size();
    auto& w=result.weights;w.fill(1.0/mixture_components);
    double old=-std::numeric_limits<double>::infinity();
    for(;result.iterations<1000;++result.iterations) {
        Values count{},logw{};for(int j=0;j<mixture_components;++j)logw[j]=fqdup::fdlibm::log(w[j]);
        double ll=0;
        for(const auto& p:pairs) {
            Values a;for(int j=0;j<mixture_components;++j)a[j]=p.log_likelihood[j]+logw[j];
            const double z=logsum(a);ll+=z;
            for(int j=0;j<mixture_components;++j)count[j]+=fqdup::fdlibm::exp(a[j]-z);
        }
        Values next;for(int j=0;j<mixture_components;++j)next[j]=std::max(count[j]/pairs.size(),1e-300);
        const double mass=sum(next);for(double& x:next)x/=mass;
        if(ll+1e-7<old)throw std::runtime_error("EM likelihood decreased");
        double delta=0;for(int j=0;j<mixture_components;++j)delta=std::max(delta,std::abs(w[j]-next[j]));
        w=next;
        if(std::isfinite(old)&&ll-old<=1e-8*(1+std::abs(old))&&delta<=1e-6){result.converged=true;++result.iterations;break;}
        old=ll;
    }
    std::vector<Values> densities;densities.reserve(pairs.size());
    for(const auto& p:pairs) { Values a;const double m=*std::max_element(p.log_likelihood.begin(),p.log_likelihood.end());for(int j=0;j<mixture_components;++j)a[j]=fqdup::fdlibm::exp(p.log_likelihood[j]-m);densities.push_back(a); }
    result.kkt_gap=std::numeric_limits<double>::infinity();
    for(;result.polish_iterations<10000;++result.polish_iterations) {
        std::vector<double> den;den.reserve(densities.size());
        for(const auto& a:densities){double d=0;for(int j=0;j<mixture_components;++j)d+=a[j]*w[j];den.push_back(d);}
        Values gradient{};
        for(std::size_t i=0;i<den.size();++i){if(!(den[i]>0))throw std::runtime_error("density underflow");for(int j=0;j<mixture_components;++j)gradient[j]+=densities[i][j]/den[i];}
        const int hi=static_cast<int>(std::max_element(gradient.begin(),gradient.end())-gradient.begin());
        result.kkt_gap=gradient[hi]-pairs.size();
        if(result.kkt_gap<=1e-6*pairs.size()){result.converged=true;break;}
        int lo=-1;for(int j=0;j<mixture_components;++j)if(w[j]>1e-14&&(lo<0||gradient[j]<gradient[lo]))lo=j;
        if(lo<0||lo==hi)break;
        std::vector<double> b;b.reserve(den.size());for(std::size_t i=0;i<den.size();++i)b.push_back((densities[i][hi]-densities[i][lo])/den[i]);
        const auto derivative=[&](double alpha){double d=0;for(double v:b)d+=v/(1+alpha*v);return d;};
        double alpha=w[lo];
        if(derivative(alpha)<0){double lower=0,upper=alpha;for(int z=0;z<55;++z){const double mid=(lower+upper)/2;if(derivative(mid)>0)lower=mid;else upper=mid;}alpha=(lower+upper)/2;}
        if(!(alpha>0))break;
        w[hi]+=alpha;w[lo]-=alpha;if(w[lo]<1e-15)w[lo]=0;
        const double mass=sum(w);for(double& x:w)x/=mass;
    }
    for(const auto& p:pairs){Values a;for(int j=0;j<mixture_components;++j)a[j]=p.log_likelihood[j]+fqdup::fdlibm::log(w[j]);result.log_likelihood+=logsum(a);}
    return result;
}
std::string prior_text(const std::array<std::optional<MixtureFit>,2>& fits) {
    std::string out="FQDUP_PAIR_MIXTURE_V1 64\n";
    for(int f=0;f<2;++f){const auto& fit=fits[f];if(fit&&!fit->converged)throw std::runtime_error("unconverged mixture");out+=std::to_string(f)+" "+std::to_string(fit?fit->training_pairs:0);for(int j=0;j<mixture_components;++j)out+=" "+shortest_double(fit?fit->weights[j]:0);out+='\n';}
    return out;
}
}
