#pragma once
#include "mixture_fit.hpp"

namespace fqdup::r8qkeep::detail {
using CountsRow=std::array<double,mixture_components>;
// Independent products/divisions use a multiple-of-four body and exact tail.
// The denominator remains the original left-to-right scalar reduction.
inline double counts_products(const CountsRow& weights,const CountsRow& density,CountsRow& out) {
    constexpr int body=(mixture_components/4)*4;
    for(int j=0;j<body;++j)out[j]=weights[j]*density[j];
    for(int j=body;j<mixture_components;++j)out[j]=weights[j]*density[j];
    double sum=0;
    for(int j=0;j<mixture_components;++j)sum+=out[j];
    return sum;
}
inline void counts_normalize(CountsRow& out,double sum) {
    constexpr int body=(mixture_components/4)*4;
    for(int j=0;j<body;++j)out[j]/=sum;
    for(int j=body;j<mixture_components;++j)out[j]/=sum;
}
}
