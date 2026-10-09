#pragma once
#include "mixture_fit.hpp"
#include "fdlibm.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace fqdup::r8qkeep::detail {
// Counts only: maximize mean log(P*w) on the simplex. Every scalar reduction
// traverses the canonical rows in order, independent of the worker count.
// A global Frank-Wolfe/KKT gap <= 1e-7 certifies the concave objective. Zero
// weights are bounds, not permanently discarded components: all 66 gradients
// are scanned at every iteration. There is no optimization iteration cap.
template<class Workers>
MixtureFit solve_counts(const std::vector<std::array<double,mixture_components>>& p,
                       const std::vector<double>& maxima,
                       const std::vector<std::size_t>& multiplicity,
                       std::size_t pairs,Workers& workers) {
    using Row=std::array<double,mixture_components>;
    constexpr double tolerance=1e-7;
    const std::size_t n=p.size();
    MixtureFit result;result.training_pairs=pairs;result.unique_likelihoods=n;
    auto& w=result.weights;w.fill(1.0/mixture_components);
    std::vector<double> den(n),candidate_den(n);
    const auto support=[](const Row& weights){
        std::vector<int> ids;for(int j=0;j<mixture_components;++j)if(weights[j]>0)ids.push_back(j);
        return ids;
    };
    const auto evaluate=[&](const Row& weights,std::vector<double>& sums){
        const auto ids=support(weights);
        workers.run(n,[&](std::size_t i,int){double d=0;for(int j:ids)d+=weights[j]*p[i][j];sums[i]=d;});
        double value=0;
        for(std::size_t i=0;i<n;++i){
            if(!(sums[i]>0)||!std::isfinite(sums[i]))return -std::numeric_limits<double>::infinity();
            value+=double(multiplicity[i])*fqdup::fdlibm::log(sums[i]);
        }
        return value/double(pairs);
    };
    const auto gradient=[&]{
        Row g{};
        workers.run((mixture_components+7)/8,[&](std::size_t block,int){
            const int first=int(block)*8,last=std::min(first+8,mixture_components);
            std::array<double,8> values{};
            for(std::size_t i=0;i<n;++i){
                const double factor=double(multiplicity[i])/den[i];
                for(int j=first;j<last;++j)values[j-first]+=factor*p[i][j];
            }
            for(int j=first;j<last;++j)g[j]=values[j-first]/double(pairs);
        },1);
        return g;
    };
    const auto normalize=[](Row& weights){double mass=0;for(double x:weights)mass+=x;for(double& x:weights)x/=mass;};
    const auto em_begin=std::chrono::steady_clock::now();
    double objective=evaluate(w,den);
    if(!std::isfinite(objective))throw std::runtime_error("counts initial density underflow");
    // A fixed warm start is not a termination cap: Newton/KKT decides convergence.
    for(int warm=0;warm<8;++warm){
        const Row g=gradient();for(int j=0;j<mixture_components;++j)w[j]*=g[j];
        normalize(w);objective=evaluate(w,den);++result.iterations;
    }
    // Start Newton on the four largest warm-start components when they cover
    // all rows. Re-entry is driven by the global KKT scan, not this initial set.
    std::array<int,mixture_components> order{};std::iota(order.begin(),order.end(),0);
    std::stable_sort(order.begin(),order.end(),[&](int a,int b){return w[a]>w[b];});
    Row sparse{};for(int k=0;k<4;++k)sparse[order[k]]=w[order[k]];
    normalize(sparse);
    const double sparse_objective=evaluate(sparse,candidate_den);
    const Row warm_gradient=gradient();
    if(*std::max_element(warm_gradient.begin(),warm_gradient.end())-1>tolerance &&
       std::isfinite(sparse_objective) && *std::min_element(candidate_den.begin(),candidate_den.end())>1e-100){w=sparse;objective=sparse_objective;den.swap(candidate_den);}
    result.em_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-em_begin).count();
    const auto newton_begin=std::chrono::steady_clock::now();
    for(;;){
        const Row g=gradient();
        const int hi=int(std::max_element(g.begin(),g.end())-g.begin());
        result.kkt_gap=(g[hi]-1)*double(pairs);
        if(g[hi]-1<=tolerance){result.converged=true;break;}
        auto ids=support(w);if(w[hi]==0){ids.push_back(hi);++result.newton_reentries;}
        std::sort(ids.begin(),ids.end());
        const int anchor=int(std::max_element(w.begin(),w.end())-w.begin());
        ids.erase(std::find(ids.begin(),ids.end(),anchor));
        Row direction{};
        // Bound-active variables with outward Newton directions stay fixed for
        // this solve. They are eligible for the next global re-entry scan.
        for(;;){
            const std::size_t m=ids.size();
            if(!m)break;
            std::vector<double> diff(m*n),h(m*m),rhs(m);
            workers.run(m,[&](std::size_t a,int){
                const int j=ids[a];rhs[a]=g[j]-g[anchor];
                for(std::size_t i=0;i<n;++i)diff[a*n+i]=(p[i][j]-p[i][anchor])/den[i];
            },1);
            workers.run(m,[&](std::size_t a,int){
                for(std::size_t b=0;b<=a;++b){
                    double value=0;
                    for(std::size_t i=0;i<n;++i)value+=double(multiplicity[i])*diff[a*n+i]*diff[b*n+i];
                    h[a*m+b]=h[b*m+a]=value/double(pairs);
                }
            },1);
            double scale=0;for(std::size_t a=0;a<m;++a)scale=std::max(scale,h[a*m+a]);
            // Damping handles coincident / nearly collinear mixture components.
            // It changes the search direction only; convergence uses the original
            // unregularized objective and gradients over every component.
            long double ridge=std::max(1e-14,scale*1e-12);
            std::vector<long double> l(m*m),step(m);
            for(;;){
                bool positive=true;std::fill(l.begin(),l.end(),0);
                for(std::size_t a=0;a<m&&positive;++a)for(std::size_t b=0;b<=a;++b){
                    long double v=h[a*m+b];if(a==b)v+=ridge;
                    for(std::size_t k=0;k<b;++k)v-=l[a*m+k]*l[b*m+k];
                    if(a==b){if(!(v>0)||!std::isfinite(v)){positive=false;break;}l[a*m+b]=std::sqrt(v);}
                    else l[a*m+b]=v/l[b*m+b];
                }
                if(positive)break;
                ridge*=100;if(!std::isfinite(ridge))throw std::runtime_error("counts Newton matrix is not finite");
            }
            for(std::size_t a=0;a<m;++a){long double v=rhs[a];for(std::size_t b=0;b<a;++b)v-=l[a*m+b]*step[b];step[a]=v/l[a*m+a];}
            for(std::size_t a=m;a-->0;){long double v=step[a];for(std::size_t b=a+1;b<m;++b)v-=l[b*m+a]*step[b];step[a]=v/l[a*m+a];}
            direction.fill(0);int outward=-1;
            for(std::size_t a=0;a<m;++a){
                const int j=ids[a];direction[j]=double(step[a]);direction[anchor]-=direction[j];
                if(w[j]==0&&direction[j]<0&&outward<0)outward=int(a);
            }
            if(outward<0)break;
            ids.erase(ids.begin()+outward);
        }
        const auto try_direction=[&](const Row& d){
            double slope=0,alpha=1;int bound=-1;
            for(int j=0;j<mixture_components;++j){
                slope+=g[j]*d[j];
                if(d[j]<0&&-w[j]/d[j]<=alpha){alpha=-w[j]/d[j];bound=j;}
            }
            if(!(slope>0)||!(alpha>0)||!std::isfinite(slope))return false;
            const double boundary=alpha;
            for(;;){
                Row next;for(int j=0;j<mixture_components;++j)next[j]=std::max(0.0,w[j]+alpha*d[j]);
                if(bound>=0&&alpha==boundary)next[bound]=0;
                normalize(next);if(next==w)return false;
                const double value=evaluate(next,candidate_den);
                if(std::isfinite(value)&&value>=objective+1e-4*alpha*slope-2e-14*(1+std::abs(objective))){
                    w=next;objective=value;den.swap(candidate_den);return true;
                }
                const double smaller=alpha*.5;if(!(smaller>0)||smaller==alpha)return false;alpha=smaller;
            }
        };
        if(!try_direction(direction)){
            // Exact one-dimensional pair exchange is a robust fallback for a
            // degenerate reduced Hessian or a Newton direction blocked by bounds.
            int lo=-1;for(int j=0;j<mixture_components;++j)if(w[j]>0&&(lo<0||g[j]<g[lo]))lo=j;
            if(lo<0||lo==hi)throw std::runtime_error("counts KKT stalled without a feasible direction");
            std::vector<double> b(n);for(std::size_t i=0;i<n;++i)b[i]=(p[i][hi]-p[i][lo])/den[i];
            const auto derivative=[&](double alpha){
                double value=0;for(std::size_t i=0;i<n;++i){const double d=1+alpha*b[i];if(!(d>0))return -std::numeric_limits<double>::infinity();value+=double(multiplicity[i])*b[i]/d;}
                return value;
            };
            double alpha=w[lo];
            if(derivative(alpha)<0){
                double lower=0,upper=alpha;
                for(;;){const double mid=lower+(upper-lower)*.5;if(mid==lower||mid==upper)break;if(derivative(mid)>0)lower=mid;else upper=mid;}
                alpha=lower+(upper-lower)*.5;
            }
            direction.fill(0);direction[hi]=alpha;direction[lo]=-alpha;
            if(!try_direction(direction))throw std::runtime_error("counts solver stalled above KKT tolerance");
        }
        ++result.polish_iterations;
    }
    result.polish_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-newton_begin).count();
    for(std::size_t i=0;i<n;++i)result.log_likelihood+=double(multiplicity[i])*(maxima[i]+fqdup::fdlibm::log(den[i]));
    return result;
}
}
