#include "mixture_fit.hpp"
#include "shortest_double.hpp"
#include "fdlibm.hpp"
#include "exact_exp_cache.hpp"
#include "worker_pool.hpp"
#include <chrono>
#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <cstring>
#include <future>
namespace fqdup::r8qkeep {
namespace {
using Values = std::array<double,mixture_components>;
struct BitsHash {std::size_t operator()(const Values* p)const {std::uint64_t hash=14695981039346656037ULL;for(double x:*p){std::uint64_t bits;std::memcpy(&bits,&x,sizeof(bits));hash^=bits;hash*=1099511628211ULL;}return hash;}};
struct BitsEqual {bool operator()(const Values* a,const Values* b)const{return std::memcmp(a->data(),b->data(),sizeof(double)*mixture_components)==0;}};
class FitWorkers {
    int threads_;
    std::unique_ptr<WorkerPool> pool_;
    std::vector<std::future<void>> jobs_;
public:
    explicit FitWorkers(int threads):threads_(threads){
        if(threads>1)pool_=std::make_unique<WorkerPool>(threads);
        jobs_.reserve(threads);
    }
    template<class Function> void run(std::size_t count,Function function){
        const int workers=std::max(1,std::min<int>(threads_,count));
        if(workers==1){for(std::size_t i=0;i<count;++i)function(i,0);return;}
        jobs_.clear();
        try {
            for(int w=0;w<workers;++w)jobs_.push_back(pool_->submit([&,w]{
                for(std::size_t i=count*w/workers;i<count*(w+1)/workers;++i)function(i,w);
            }));
            for(auto& job:jobs_)job.get();
        }catch(...){
            for(auto& job:jobs_)if(job.valid())job.wait();
            throw;
        }
    }
};
double sum(const Values& values) { double out=0;for(double v:values)out+=v;return out; }
double logsum(const Values& values, ExactExpCache& exp) {
    const double m=*std::max_element(values.begin(),values.end());
    if(!std::isfinite(m))throw std::runtime_error("all mixture components impossible");
    double s=0;for(double x:values)s+=exp(x-m);
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
LikelihoodAccumulator::LikelihoodAccumulator(){
    const auto likelihood=[](double d,bool agree){return fqdup::fdlibm::log(std::max(1e-15,std::min(1-1e-15,agree?1-d:d)));};
    for(int j=1;j<mixture_components;++j){const double r=double(j-1)/64,e=.75*(r*r),d=2*e-4*e*e/3;flat_[j]={likelihood(d,true),likelihood(d,false)};}
}
void LikelihoodAccumulator::add(const std::string& id,int fold,double e1,double e2,bool agree){
    auto [it,inserted]=slots_.emplace(id,pairs_.size());if(inserted){PairLikelihood p;p.id=id;p.fold=fold;pairs_.push_back(std::move(p));}
    auto& p=pairs_[it->second];if(p.fold!=fold)throw std::runtime_error("fold mismatch");
    const double d=e1+e2-4*e1*e2/3;
    p.log_likelihood[0]+=fqdup::fdlibm::log(std::max(1e-15,std::min(1-1e-15,agree?1-d:d)));
    for(int j=1;j<mixture_components;++j)p.log_likelihood[j]+=flat_[j][agree?0:1];
}
std::optional<MixtureFit> fit_mixture(const std::vector<PairLikelihood>& pairs,int threads) {
    std::vector<const PairLikelihood*> view;view.reserve(pairs.size());for(const auto& pair:pairs)view.push_back(&pair);
    return fit_mixture_views(view,threads);
}
template<bool Counts>
static std::optional<MixtureFit> fit_mixture_impl(const std::vector<const PairLikelihood*>& pairs,int threads) {
    if(pairs.empty())return std::nullopt;
    threads=std::max(1,threads);
    std::unordered_map<const Values*,std::size_t,BitsHash,BitsEqual> slots;
    std::vector<const Values*> keys;std::vector<std::size_t> group;group.reserve(pairs.size());
    for(const auto* pair:pairs){const auto* key=&pair->log_likelihood;auto [it,inserted]=slots.emplace(key,keys.size());if(inserted)keys.push_back(key);group.push_back(it->second);}
    std::vector<std::size_t> multiplicity;
    if constexpr(Counts){
        // Numeric ordering of the exact IEEE bit patterns, independent of IDs,
        // arrival order, hash iteration, and worker scheduling.
        std::sort(keys.begin(),keys.end(),[](const Values* a,const Values* b){
            for(int j=0;j<mixture_components;++j){std::uint64_t x,y;std::memcpy(&x,&(*a)[j],sizeof(x));std::memcpy(&y,&(*b)[j],sizeof(y));if(x!=y)return x<y;}
            return false;
        });
        slots.clear();
        for(std::size_t i=0;i<keys.size();++i)slots.emplace(keys[i],i);
        multiplicity.resize(keys.size());
        for(const auto* pair:pairs){
            if(pair->fold!=pairs.front()->fold)throw std::invalid_argument("counts fit requires one training fold");
            ++multiplicity[slots.at(&pair->log_likelihood)];
        }
        group.resize(keys.size());for(std::size_t i=0;i<group.size();++i)group[i]=i;
    }
    // Tiny fits stay serial; term arithmetic and ordered reductions are identical.
    if(keys.size()<64)threads=1;
    FitWorkers workers(threads);
    std::vector<ExactExpCache> exp(threads);
    struct Terms {Values responsibility;double log_sum=0;};
    std::vector<Terms> terms(keys.size());
    MixtureFit result;result.training_pairs=pairs.size();
    result.unique_likelihoods=keys.size();
    auto& w=result.weights;w.fill(1.0/mixture_components);
    // Counts reuses these row-normalized likelihoods for EM and polishing.
    // Full deliberately keeps its original log-domain arithmetic and ordering.
    std::vector<Values> densities(keys.size());
    std::vector<double> row_max;
    if constexpr(Counts){
        row_max.resize(keys.size());
        workers.run(keys.size(),[&](std::size_t i,int worker){
            const auto& p=*keys[i];
            const double m=*std::max_element(p.begin(),p.end());
            if(!std::isfinite(m))throw std::runtime_error("all mixture components impossible");
            row_max[i]=m;
            for(int j=0;j<mixture_components;++j)densities[i][j]=exp[worker](p[j]-m);
        });
    }
    const auto em_begin=std::chrono::steady_clock::now();
    double old=-std::numeric_limits<double>::infinity();
    for(;result.iterations<1000;++result.iterations) {
        Values count{},logw{};for(int j=0;j<mixture_components;++j)logw[j]=fqdup::fdlibm::log(w[j]);
        double ll=0;
        workers.run(keys.size(),[&](std::size_t i,int worker){
            auto& t=terms[i];
            if constexpr(Counts){
                double s=0;
                for(int j=0;j<mixture_components;++j){t.responsibility[j]=w[j]*densities[i][j];s+=t.responsibility[j];}
                if(s>0){
                    t.log_sum=row_max[i]+fqdup::fdlibm::log(s);
                    for(double& r:t.responsibility)r/=s;
                    return;
                }
                // A zero scaled denominator must use the original log path.
            }
            Values a;for(int j=0;j<mixture_components;++j)a[j]=(*keys[i])[j]+logw[j];
            t.log_sum=logsum(a,exp[worker]);
            for(int j=0;j<mixture_components;++j)t.responsibility[j]=exp[worker](a[j]-t.log_sum);
        });
        // Every floating-point reduction still visits every original pair in order.
        for(std::size_t index:group){const auto& t=terms[index];if constexpr(Counts){const double n=multiplicity[index];ll+=n*t.log_sum;for(int j=0;j<mixture_components;++j)count[j]+=n*t.responsibility[j];}else{ll+=t.log_sum;for(int j=0;j<mixture_components;++j)count[j]+=t.responsibility[j];}}
        Values next;for(int j=0;j<mixture_components;++j)next[j]=std::max(count[j]/pairs.size(),1e-300);
        const double mass=sum(next);for(double& x:next)x/=mass;
        if(ll+1e-7<old)throw std::runtime_error("EM likelihood decreased");
        double delta=0;for(int j=0;j<mixture_components;++j)delta=std::max(delta,std::abs(w[j]-next[j]));
        w=next;
        if(std::isfinite(old)&&ll-old<=1e-8*(1+std::abs(old))&&delta<=1e-6){result.converged=true;++result.iterations;break;}
        old=ll;
    }
    result.em_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-em_begin).count();
    if constexpr(!Counts)workers.run(keys.size(),[&](std::size_t i,int worker){const auto& p=*keys[i];const double m=*std::max_element(p.begin(),p.end());for(int j=0;j<mixture_components;++j)densities[i][j]=exp[worker](p[j]-m);});
    result.kkt_gap=std::numeric_limits<double>::infinity();
    std::vector<double> table_den(densities.size()),den(group.size()),b(group.size());
    const auto polish_begin=std::chrono::steady_clock::now();
    for(;result.polish_iterations<10000;++result.polish_iterations) {
        workers.run(densities.size(),[&](std::size_t i,int){double d=0;for(int j=0;j<mixture_components;++j)d+=densities[i][j]*w[j];table_den[i]=d;});
        for(std::size_t i=0;i<group.size();++i)den[i]=table_den[group[i]];
        Values gradient{};
        for(std::size_t i=0;i<den.size();++i){if(!(den[i]>0))throw std::runtime_error("density underflow");const auto& density=densities[group[i]];for(int j=0;j<mixture_components;++j){if constexpr(Counts)gradient[j]+=multiplicity[i]*(density[j]/den[i]);else gradient[j]+=density[j]/den[i];}}
        const int hi=static_cast<int>(std::max_element(gradient.begin(),gradient.end())-gradient.begin());
        result.kkt_gap=gradient[hi]-pairs.size();
        if(result.kkt_gap<=1e-6*pairs.size()){result.converged=true;break;}
        int lo=-1;for(int j=0;j<mixture_components;++j)if(w[j]>1e-14&&(lo<0||gradient[j]<gradient[lo]))lo=j;
        if(lo<0||lo==hi)break;
        for(std::size_t i=0;i<den.size();++i)b[i]=(densities[group[i]][hi]-densities[group[i]][lo])/den[i];
        const auto derivative=[&](double alpha){double d=0;if constexpr(Counts){for(std::size_t i=0;i<b.size();++i){const double v=b[i];d+=multiplicity[i]*(v/(1+alpha*v));}}else{for(double v:b)d+=v/(1+alpha*v);}return d;};
        double alpha=w[lo];
        if(derivative(alpha)<0){double lower=0,upper=alpha;for(int z=0;z<55;++z){const double mid=(lower+upper)/2;if(derivative(mid)>0)lower=mid;else upper=mid;}alpha=(lower+upper)/2;}
        if(!(alpha>0))break;
        w[hi]+=alpha;w[lo]-=alpha;if(w[lo]<1e-15)w[lo]=0;
        const double mass=sum(w);for(double& x:w)x/=mass;
    }
    result.polish_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-polish_begin).count();
    workers.run(keys.size(),[&](std::size_t i,int worker){
        if constexpr(Counts){
            double s=0;for(int j=0;j<mixture_components;++j)s+=w[j]*densities[i][j];
            if(s>0){terms[i].log_sum=row_max[i]+fqdup::fdlibm::log(s);return;}
        }
        Values a;for(int j=0;j<mixture_components;++j)a[j]=(*keys[i])[j]+fqdup::fdlibm::log(w[j]);
        terms[i].log_sum=logsum(a,exp[worker]);
    });
    for(auto index:group){if constexpr(Counts)result.log_likelihood+=multiplicity[index]*terms[index].log_sum;else result.log_likelihood+=terms[index].log_sum;}
    for(const auto& cache:exp){result.exp_calls+=cache.calls();result.exp_cache_hits+=cache.hits();}
    return result;
}
std::optional<MixtureFit> fit_mixture_views(const std::vector<const PairLikelihood*>& pairs,int threads){return fit_mixture_impl<false>(pairs,threads);}
std::optional<MixtureFit> fit_mixture_counts(const std::vector<const PairLikelihood*>& pairs,int threads){return fit_mixture_impl<true>(pairs,threads);}
std::string prior_text(const std::array<std::optional<MixtureFit>,2>& fits) {
    std::string out="FQDUP_PAIR_MIXTURE_V1 64\n";
    for(int f=0;f<2;++f){const auto& fit=fits[f];if(fit&&!fit->converged)throw std::runtime_error("unconverged mixture");out+=std::to_string(f)+" "+std::to_string(fit?fit->training_pairs:0);for(int j=0;j<mixture_components;++j)out+=" "+shortest_double(fit?fit->weights[j]:0);out+='\n';}
    return out;
}
}
