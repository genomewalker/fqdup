#include "../src/merge_r8qkeep/mixture_fit.hpp"
#include "../src/merge_r8qkeep/worker_pool.hpp"
#include "../src/merge_r8qkeep/counts_row.hpp"
#include "../src/merge_r8qkeep/fdlibm.hpp"
#include <cfenv>
#include <random>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <fstream>
#include <future>
#include <chrono>
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
    require(std::memcmp(a.weights.data(),b.weights.data(),sizeof(a.weights))==0 &&
        std::memcmp(&a.kkt_gap,&b.kkt_gap,sizeof(double))==0 &&
        std::memcmp(&a.log_likelihood,&b.log_likelihood,sizeof(double))==0,"fit bit-pattern dependence");
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
static int benchmark_workers(const char* census_path){
    std::ifstream input(census_path);
    require(bool(input),"cannot open worker benchmark census");
    const auto pairs=read_census_likelihood(input,true);
    std::array<std::vector<const PairLikelihood*>,2> training;
    for(const auto& p:pairs)training[p.fold].push_back(&p);
    std::array<std::optional<MixtureFit>,2> baseline;
    for(int workers:{1,2,4,4,2,1}){
        const auto begin=std::chrono::steady_clock::now();
        std::array<std::future<std::optional<MixtureFit>>,2> jobs;
        for(int fold=0;fold<2;++fold)
            jobs[fold]=std::async(std::launch::async,[&,fold,workers]{
                return fit_mixture_counts(training[fold],workers);
            });
        std::array<std::optional<MixtureFit>,2> fits;
        for(int fold=0;fold<2;++fold)fits[fold]=jobs[fold].get();
        const double wall=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
        for(int fold=0;fold<2;++fold){
            require(bool(fits[fold]),"empty worker benchmark fit");
            if(!baseline[fold])baseline[fold]=fits[fold];
            identical(*baseline[fold],*fits[fold]);
        }
        require(prior_text(baseline)==prior_text(fits),"worker benchmark prior differs");
        std::cout<<"FIT_WORKERS workers_per_fit="<<workers<<" wall="<<wall
                 <<" keys="<<fits[0]->unique_likelihoods<<","<<fits[1]->unique_likelihoods
                 <<" iterations="<<fits[0]->iterations<<","<<fits[1]->iterations
                 <<" exact=1"<<std::endl;
    }
    for(int repeat=0;repeat<2;++repeat){
        const auto begin=std::chrono::steady_clock::now();
        std::array<std::optional<MixtureFit>,2> fits;
        {
            WorkerPool pool(6);
            std::array<std::future<std::optional<MixtureFit>>,2> jobs;
            for(int fold=0;fold<2;++fold)jobs[fold]=std::async(std::launch::async,[&,fold]{
                return fit_mixture_counts(training[fold],7,&pool);
            });
            for(int fold=0;fold<2;++fold)fits[fold]=jobs[fold].get();
        }
        const double wall=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
        for(int fold=0;fold<2;++fold){require(bool(fits[fold]),"empty shared fit");identical(*baseline[fold],*fits[fold]);}
        require(prior_text(baseline)==prior_text(fits),"shared worker prior differs");
        std::cout<<"FIT_SHARED total_workers=8 logical_workers_per_fit=7 wall="<<wall<<" exact=1"<<std::endl;
    }
    return 0;
}
static void row_terms_parity(){
    const int original_round=std::fegetround();
    std::mt19937_64 rng(2530);
    for(int mode:{FE_TONEAREST,FE_DOWNWARD,FE_UPWARD,FE_TOWARDZERO}){
        require(std::fesetround(mode)==0,"set row-oracle rounding mode");
        for(int trial=0;trial<256;++trial){
            detail::CountsRow weights,density,expected,actual;
            for(int j=0;j<mixture_components;++j){
                weights[j]=double(1+rng()%1000)/1001.0;
                if(trial%3==0)weights[j]=std::ldexp(weights[j],-int(rng()%1000));
                switch((j+trial)%8){
                    case 0:density[j]=0;break;
                    case 1:density[j]=std::numeric_limits<double>::denorm_min()*double(1+rng()%1000);break;
                    case 2:density[j]=std::numeric_limits<double>::min();break;
                    case 3:density[j]=1;break;
                    case 4:density[j]=std::ldexp(double(1+rng()%1000)/1001.0,-int(rng()%1000));break;
                    case 5:density[j]=std::nextafter(1.0,0.0);break;
                    case 6:density[j]=std::nextafter(std::numeric_limits<double>::min(),0.0);break;
                    default:density[j]=.5;
                }
                if(trial==0)density[j]=0;
            }
            // Literal previous interleaved loop, including every ordered add.
            double expected_sum=0;
            for(int j=0;j<mixture_components;++j){
                expected[j]=weights[j]*density[j];expected_sum+=expected[j];
            }
            const double actual_sum=detail::counts_products(weights,density,actual);
            require(std::memcmp(expected.data(),actual.data(),sizeof(expected))==0 &&
                    std::memcmp(&expected_sum,&actual_sum,sizeof(double))==0,"counts products/sum bits");
            if(expected_sum>0){
                for(double& value:expected)value/=expected_sum;
                detail::counts_normalize(actual,actual_sum);
                require(std::memcmp(expected.data(),actual.data(),sizeof(expected))==0,"counts division bits");
            }
        }
    }
    require(std::fesetround(original_round)==0,"restore row-oracle rounding mode");
}
int main(int argc,char** argv){
    if(argc==2)return benchmark_workers(argv[1]);
    require(argc==1,"unexpected counts test arguments");
    row_terms_parity();
    // A partial final 64-row chunk exercises bounded dynamic scheduling and unwind.
    std::vector<PairLikelihood> parallel(513);
    for(std::size_t i=0;i<parallel.size();++i){
        parallel[i].id=std::to_string(i);
        for(int j=0;j<mixture_components;++j)
            parallel[i].log_likelihood[j]=(j==int(i%3)?0.0:-20.0)-.001*i;
    }
    const auto parallel_view=view(parallel);
    const auto serial_counts=fit_mixture_counts(parallel_view,1),serial_full=fit_mixture_views(parallel_view,1);
    for(int threads:{2,4,8}){
        identical(*serial_counts,*fit_mixture_counts(parallel_view,threads));
        identical(*serial_full,*fit_mixture_views(parallel_view,threads));
    }
    for(bool counts:{false,true})for(std::size_t bad:{std::size_t(0),parallel.size()-1}){
        auto poison=parallel;poison[bad].log_likelihood.fill(-std::numeric_limits<double>::infinity());
        bool worker_error=false;
        try{if(counts)fit_mixture_counts(view(poison),8);else fit_mixture_views(view(poison),8);}
        catch(const std::runtime_error&){worker_error=true;}
        require(worker_error,"caller or worker exception was lost");
    }
    identical(*fit_mixture_counts(parallel_view,1),*fit_mixture_counts(parallel_view,8));
    // One caller can fail while another continues using the borrowed pool.
    {
        WorkerPool pool(2);
        auto good=std::async(std::launch::async,[&]{return fit_mixture_counts(parallel_view,3,&pool);});
        auto bad=std::async(std::launch::async,[&]{
            auto poison=parallel;poison.back().log_likelihood.fill(-std::numeric_limits<double>::infinity());
            return fit_mixture_counts(view(poison),3,&pool);
        });
        bool threw=false;try{bad.get();}catch(const std::runtime_error&){threw=true;}
        require(threw,"shared fit exception was lost");identical(*serial_counts,*good.get());
        identical(*serial_counts,*fit_mixture_counts(parallel_view,3,&pool));
    }
    // Normal, subnormal and zero densities exercise the reordered independent
    // accumulator schedule against the literal serial reduction path.
    std::vector<PairLikelihood> deep(257);
    for(std::size_t i=0;i<deep.size();++i){deep[i].id=std::to_string(i);
        for(int j=0;j<mixture_components;++j)
            deep[i].log_likelihood[j]=(j==int(i%4)?0.0:-20.0-100.0*(j%9))-.001*i;
    }
    const auto deep_view=view(deep);
    const auto deep_serial=fit_mixture_counts(deep_view,1);
    for(int threads:{2,4,8})identical(*deep_serial,*fit_mixture_counts(deep_view,threads));
    {
        WorkerPool pool(6);
        auto first=std::async(std::launch::async,[&]{return fit_mixture_counts(deep_view,7,&pool);});
        auto second=std::async(std::launch::async,[&]{return fit_mixture_counts(deep_view,7,&pool);});
        identical(*deep_serial,*first.get());identical(*deep_serial,*second.get());
    }
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
    // Components outside the four-component initial support must re-enter.
    std::vector<PairLikelihood> reentry;
    for(int kind=0;kind<6;++kind)for(int copy=0;copy<10*(kind+1);++copy){
        PairLikelihood p;p.id=std::to_string(reentry.size());p.log_likelihood.fill(-40);
        for(int j=0;j<6;++j)p.log_likelihood[j]=-3;
        p.log_likelihood[kind]=0;reentry.push_back(p);
    }
    const auto entered=fit_mixture_counts(view(reentry));
    require(entered->newton_reentries>0,"inactive components never re-entered");
    repeated_objective(reentry,*entered);
    for(int j=0;j<6;++j){const double a=std::exp(-3.0),expected=(double(j+1)/21*(1+5*a)-a)/(1-a);
        require(std::abs(entered->weights[j]-expected)<1e-6,"re-entry optimum");}
    std::reverse(reentry.begin(),reentry.end());
    identical(*entered,*fit_mixture_counts(view(reentry),8));

    // Per-pair statistics preserve duplicate IDs/folds and nominal addition order;
    // flat components use A*log(1-d)+D*log(d), not pooled base counts.
    LikelihoodAccumulator stats(true),expanded_stats;
    for(int k=0;k<120;++k)for(int id=0;id<6;++id){
        const bool agree=k%(id+3)!=0;
        const std::string name=std::to_string(id);
        stats.add(name,1,.01,.02,agree);expanded_stats.add(name,1,.01,.02,agree);
    }
    const auto counted=stats.take();auto manual=expanded_stats.take();
    require(counted.size()==6,"duplicate ID statistics split");
    for(std::size_t i=0;i<counted.size();++i){
        const auto& p=counted[i];
        require(p.sufficient_statistics&&p.agreements+p.disagreements==120,"pair A/D counts");
        require(p.log_likelihood[0]==manual[i].log_likelihood[0],"nominal accumulation changed");
        for(int j=1;j<mixture_components;++j){
            const double r=double(j-1)/64,e=.75*(r*r),d=2*e-4*e*e/3;
            manual[i].log_likelihood[j]=double(p.agreements)*fqdup::fdlibm::log(std::max(1e-15,std::min(1-1e-15,1-d)))+
                double(p.disagreements)*fqdup::fdlibm::log(std::max(1e-15,std::min(1-1e-15,d)));
        }
    }
    const auto stat_fit=fit_mixture_counts(view(counted));
    identical(*stat_fit,*fit_mixture_counts(view(manual),8));
    repeated_objective(manual,*stat_fit);
    LikelihoodAccumulator wrong_fold(true);wrong_fold.add("duplicate",0,.01,.02,true);
    bool fold_error=false;try{wrong_fold.add("duplicate",1,.01,.02,false);}catch(const std::runtime_error&){fold_error=true;}
    require(fold_error,"statistics accepted duplicate ID fold mismatch");
    std::cout<<"TASK25_COUNTS_PASS\n";
}
