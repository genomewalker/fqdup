#pragma once
// V8-ported fdlibm routines preserve prior-fit arithmetic for byte reproducibility.
namespace fqdup::fdlibm {
double exp(double x);
double log(double x);
double log1p(double x);
double expm1(double x);
double pow(double x,double y);
}
