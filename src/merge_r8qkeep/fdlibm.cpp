// The following is adapted from fdlibm (http://www.netlib.org/fdlibm).
//
// ====================================================
// Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
//
// Developed at SunSoft, a Sun Microsystems, Inc. business.
// Permission to use, copy, modify, and distribute this
// software is freely granted, provided that this notice
// is preserved.
// ====================================================
//
// The original source code covered by the above license above has been
// modified significantly by Google Inc.
// Copyright 2016 the V8 project authors. All rights reserved.

// Vendored from Node v22.20.0 / V8 ieee754.cc; algorithm bodies retained.
#include "fdlibm.hpp"
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <type_traits>
namespace fqdup::fdlibm {
namespace base {
template<class To,class From> To bit_cast(const From& value) { static_assert(sizeof(To)==sizeof(From));To out;std::memcpy(&out,&value,sizeof(out));return out; }
inline double Divide(double numerator,double denominator) { return numerator/denominator; }
}
/* Get two 32 bit ints from a double.  */

#define EXTRACT_WORDS(ix0, ix1, d)               \
  do {                                           \
    uint64_t bits = base::bit_cast<uint64_t>(d); \
    (ix0) = bits >> 32;                          \
    (ix1) = bits & 0xFFFFFFFFu;                  \
  } while (false)

/* Get the more significant 32 bit int from a double.  */

#define GET_HIGH_WORD(i, d)                      \
  do {                                           \
    uint64_t bits = base::bit_cast<uint64_t>(d); \
    (i) = bits >> 32;                            \
  } while (false)

/* Get the less significant 32 bit int from a double.  */

#define GET_LOW_WORD(i, d)                       \
  do {                                           \
    uint64_t bits = base::bit_cast<uint64_t>(d); \
    (i) = bits & 0xFFFFFFFFu;                    \
  } while (false)

/* Set a double from two 32 bit ints.  */

#define INSERT_WORDS(d, ix0, ix1)             \
  do {                                        \
    uint64_t bits = 0;                        \
    bits |= static_cast<uint64_t>(ix0) << 32; \
    bits |= static_cast<uint32_t>(ix1);       \
    (d) = base::bit_cast<double>(bits);       \
  } while (false)

/* Set the more significant 32 bits of a double from an int.  */

#define SET_HIGH_WORD(d, v)                      \
  do {                                           \
    uint64_t bits = base::bit_cast<uint64_t>(d); \
    bits &= 0x0000'0000'FFFF'FFFF;               \
    bits |= static_cast<uint64_t>(v) << 32;      \
    (d) = base::bit_cast<double>(bits);          \
  } while (false)

/* Set the less significant 32 bits of a double from an int.  */

#define SET_LOW_WORD(d, v)                       \
  do {                                           \
    uint64_t bits = base::bit_cast<uint64_t>(d); \
    bits &= 0xFFFF'FFFF'0000'0000;               \
    bits |= static_cast<uint32_t>(v);            \
    (d) = base::bit_cast<double>(bits);          \
  } while (false)

double exp(double x) {
  static const double
      one = 1.0,
      halF[2] = {0.5, -0.5},
      o_threshold = 7.09782712893383973096e+02,  /* 0x40862E42, 0xFEFA39EF */
      u_threshold = -7.45133219101941108420e+02, /* 0xC0874910, 0xD52D3051 */
      ln2HI[2] = {6.93147180369123816490e-01,    /* 0x3FE62E42, 0xFEE00000 */
                  -6.93147180369123816490e-01},  /* 0xBFE62E42, 0xFEE00000 */
      ln2LO[2] = {1.90821492927058770002e-10,    /* 0x3DEA39EF, 0x35793C76 */
                  -1.90821492927058770002e-10},  /* 0xBDEA39EF, 0x35793C76 */
      invln2 = 1.44269504088896338700e+00,       /* 0x3FF71547, 0x652B82FE */
      P1 = 1.66666666666666019037e-01,           /* 0x3FC55555, 0x5555553E */
      P2 = -2.77777777770155933842e-03,          /* 0xBF66C16C, 0x16BEBD93 */
      P3 = 6.61375632143793436117e-05,           /* 0x3F11566A, 0xAF25DE2C */
      P4 = -1.65339022054652515390e-06,          /* 0xBEBBBD41, 0xC5D26BF1 */
      P5 = 4.13813679705723846039e-08,           /* 0x3E663769, 0x72BEA4D0 */
      E = 2.718281828459045;                     /* 0x4005BF0A, 0x8B145769 */

  static volatile double
      huge = 1.0e+300,
      twom1000 = 9.33263618503218878990e-302, /* 2**-1000=0x01700000,0*/
      two1023 = 8.988465674311579539e307;     /* 0x1p1023 */

  double y, hi = 0.0, lo = 0.0, c, t, twopk;
  int32_t k = 0, xsb;
  uint32_t hx;

  GET_HIGH_WORD(hx, x);
  xsb = (hx >> 31) & 1; /* sign bit of x */
  hx &= 0x7FFFFFFF;     /* high word of |x| */

  /* filter out non-finite argument */
  if (hx >= 0x40862E42) { /* if |x|>=709.78... */
    if (hx >= 0x7FF00000) {
      uint32_t lx;
      GET_LOW_WORD(lx, x);
      if (((hx & 0xFFFFF) | lx) != 0)
        return x + x; /* NaN */
      else
        return (xsb == 0) ? x : 0.0; /* exp(+-inf)={inf,0} */
    }
    if (x > o_threshold) return huge * huge;         /* overflow */
    if (x < u_threshold) return twom1000 * twom1000; /* underflow */
  }

  /* argument reduction */
  if (hx > 0x3FD62E42) {   /* if  |x| > 0.5 ln2 */
    if (hx < 0x3FF0A2B2) { /* and |x| < 1.5 ln2 */
      /* TODO(rtoy): We special case exp(1) here to return the correct
       * value of E, as the computation below would get the last bit
       * wrong. We should probably fix the algorithm instead.
       */
      if (x == 1.0) return E;
      hi = x - ln2HI[xsb];
      lo = ln2LO[xsb];
      k = 1 - xsb - xsb;
    } else {
      k = static_cast<int>(invln2 * x + halF[xsb]);
      t = k;
      hi = x - t * ln2HI[0]; /* t*ln2HI is exact here */
      lo = t * ln2LO[0];
    }
    x = hi - lo;
  } else if (hx < 0x3E300000) {         /* when |x|<2**-28 */
    if (huge + x > one) return one + x; /* trigger inexact */
  } else {
    k = 0;
  }

  /* x is now in primary range */
  t = x * x;
  if (k >= -1021) {
    INSERT_WORDS(
        twopk,
        0x3FF00000 + static_cast<int32_t>(static_cast<uint32_t>(k) << 20), 0);
  } else {
    INSERT_WORDS(twopk, 0x3FF00000 + (static_cast<uint32_t>(k + 1000) << 20),
                 0);
  }
  c = x - t * (P1 + t * (P2 + t * (P3 + t * (P4 + t * P5))));
  if (k == 0) {
    return one - ((x * c) / (c - 2.0) - x);
  } else {
    y = one - ((lo - (x * c) / (2.0 - c)) - hi);
  }
  if (k >= -1021) {
    if (k == 1024) return y * 2.0 * two1023;
    return y * twopk;
  } else {
    return y * twopk * twom1000;
  }
}

double log(double x) {
  static const double                      /* -- */
      ln2_hi = 6.93147180369123816490e-01, /* 3fe62e42 fee00000 */
      ln2_lo = 1.90821492927058770002e-10, /* 3dea39ef 35793c76 */
      two54 = 1.80143985094819840000e+16,  /* 43500000 00000000 */
      Lg1 = 6.666666666666735130e-01,      /* 3FE55555 55555593 */
      Lg2 = 3.999999999940941908e-01,      /* 3FD99999 9997FA04 */
      Lg3 = 2.857142874366239149e-01,      /* 3FD24924 94229359 */
      Lg4 = 2.222219843214978396e-01,      /* 3FCC71C5 1D8E78AF */
      Lg5 = 1.818357216161805012e-01,      /* 3FC74664 96CB03DE */
      Lg6 = 1.531383769920937332e-01,      /* 3FC39A09 D078C69F */
      Lg7 = 1.479819860511658591e-01;      /* 3FC2F112 DF3E5244 */

  static const double zero = 0.0;

  double hfsq, f, s, z, R, w, t1, t2, dk;
  int32_t k, hx, i, j;
  uint32_t lx;

  EXTRACT_WORDS(hx, lx, x);

  k = 0;
  if (hx < 0x00100000) { /* x < 2**-1022  */
    if (((hx & 0x7FFFFFFF) | lx) == 0) {
      return -std::numeric_limits<double>::infinity(); /* log(+-0)=-inf */
    }
    if (hx < 0) {
      return std::numeric_limits<double>::signaling_NaN(); /* log(-#) = NaN */
    }
    k -= 54;
    x *= two54; /* subnormal number, scale up x */
    GET_HIGH_WORD(hx, x);
  }
  if (hx >= 0x7FF00000) return x + x;
  k += (hx >> 20) - 1023;
  hx &= 0x000FFFFF;
  i = (hx + 0x95F64) & 0x100000;
  SET_HIGH_WORD(x, hx | (i ^ 0x3FF00000)); /* normalize x or x/2 */
  k += (i >> 20);
  f = x - 1.0;
  if ((0x000FFFFF & (2 + hx)) < 3) { /* -2**-20 <= f < 2**-20 */
    if (f == zero) {
      if (k == 0) {
        return zero;
      } else {
        dk = static_cast<double>(k);
        return dk * ln2_hi + dk * ln2_lo;
      }
    }
    R = f * f * (0.5 - 0.33333333333333333 * f);
    if (k == 0) {
      return f - R;
    } else {
      dk = static_cast<double>(k);
      return dk * ln2_hi - ((R - dk * ln2_lo) - f);
    }
  }
  s = f / (2.0 + f);
  dk = static_cast<double>(k);
  z = s * s;
  i = hx - 0x6147A;
  w = z * z;
  j = 0x6B851 - hx;
  t1 = w * (Lg2 + w * (Lg4 + w * Lg6));
  t2 = z * (Lg1 + w * (Lg3 + w * (Lg5 + w * Lg7)));
  i |= j;
  R = t2 + t1;
  if (i > 0) {
    hfsq = 0.5 * f * f;
    if (k == 0)
      return f - (hfsq - s * (hfsq + R));
    else
      return dk * ln2_hi - ((hfsq - (s * (hfsq + R) + dk * ln2_lo)) - f);
  } else {
    if (k == 0)
      return f - s * (f - R);
    else
      return dk * ln2_hi - ((s * (f - R) - dk * ln2_lo) - f);
  }
}

double log1p(double x) {
  static const double                      /* -- */
      ln2_hi = 6.93147180369123816490e-01, /* 3fe62e42 fee00000 */
      ln2_lo = 1.90821492927058770002e-10, /* 3dea39ef 35793c76 */
      two54 = 1.80143985094819840000e+16,  /* 43500000 00000000 */
      Lp1 = 6.666666666666735130e-01,      /* 3FE55555 55555593 */
      Lp2 = 3.999999999940941908e-01,      /* 3FD99999 9997FA04 */
      Lp3 = 2.857142874366239149e-01,      /* 3FD24924 94229359 */
      Lp4 = 2.222219843214978396e-01,      /* 3FCC71C5 1D8E78AF */
      Lp5 = 1.818357216161805012e-01,      /* 3FC74664 96CB03DE */
      Lp6 = 1.531383769920937332e-01,      /* 3FC39A09 D078C69F */
      Lp7 = 1.479819860511658591e-01;      /* 3FC2F112 DF3E5244 */

  static const double zero = 0.0;

  double hfsq, f, c, s, z, R, u;
  int32_t k, hx, hu, ax;

  GET_HIGH_WORD(hx, x);
  ax = hx & 0x7FFFFFFF;

  k = 1;
  if (hx < 0x3FDA827A) {    /* 1+x < sqrt(2)+ */
    if (ax >= 0x3FF00000) { /* x <= -1.0 */
      if (x == -1.0)
        return -std::numeric_limits<double>::infinity(); /* log1p(-1)=+inf */
      else
        return std::numeric_limits<double>::signaling_NaN();  // log1p(x<-1)=NaN
    }
    if (ax < 0x3E200000) {    /* |x| < 2**-29 */
      if (two54 + x > zero    /* raise inexact */
          && ax < 0x3C900000) /* |x| < 2**-54 */
        return x;
      else
        return x - x * x * 0.5;
    }
    if (hx > 0 || hx <= static_cast<int32_t>(0xBFD2BEC4)) {
      k = 0;
      f = x;
      hu = 1;
    } /* sqrt(2)/2- <= 1+x < sqrt(2)+ */
  }
  if (hx >= 0x7FF00000) return x + x;
  if (k != 0) {
    if (hx < 0x43400000) {
      u = 1.0 + x;
      GET_HIGH_WORD(hu, u);
      k = (hu >> 20) - 1023;
      c = (k > 0) ? 1.0 - (u - x) : x - (u - 1.0); /* correction term */
      c /= u;
    } else {
      u = x;
      GET_HIGH_WORD(hu, u);
      k = (hu >> 20) - 1023;
      c = 0;
    }
    hu &= 0x000FFFFF;
    /*
     * The approximation to sqrt(2) used in thresholds is not
     * critical.  However, the ones used above must give less
     * strict bounds than the one here so that the k==0 case is
     * never reached from here, since here we have committed to
     * using the correction term but don't use it if k==0.
     */
    if (hu < 0x6A09E) {                  /* u ~< sqrt(2) */
      SET_HIGH_WORD(u, hu | 0x3FF00000); /* normalize u */
    } else {
      k += 1;
      SET_HIGH_WORD(u, hu | 0x3FE00000); /* normalize u/2 */
      hu = (0x00100000 - hu) >> 2;
    }
    f = u - 1.0;
  }
  hfsq = 0.5 * f * f;
  if (hu == 0) { /* |f| < 2**-20 */
    if (f == zero) {
      if (k == 0) {
        return zero;
      } else {
        c += k * ln2_lo;
        return k * ln2_hi + c;
      }
    }
    R = hfsq * (1.0 - 0.66666666666666666 * f);
    if (k == 0)
      return f - R;
    else
      return k * ln2_hi - ((R - (k * ln2_lo + c)) - f);
  }
  s = f / (2.0 + f);
  z = s * s;
  R = z * (Lp1 +
           z * (Lp2 + z * (Lp3 + z * (Lp4 + z * (Lp5 + z * (Lp6 + z * Lp7))))));
  if (k == 0)
    return f - (hfsq - s * (hfsq + R));
  else
    return k * ln2_hi - ((hfsq - (s * (hfsq + R) + (k * ln2_lo + c))) - f);
}

double expm1(double x) {
  static const double
      one = 1.0,
      tiny = 1.0e-300,
      o_threshold = 7.09782712893383973096e+02, /* 0x40862E42, 0xFEFA39EF */
      ln2_hi = 6.93147180369123816490e-01,      /* 0x3FE62E42, 0xFEE00000 */
      ln2_lo = 1.90821492927058770002e-10,      /* 0x3DEA39EF, 0x35793C76 */
      invln2 = 1.44269504088896338700e+00,      /* 0x3FF71547, 0x652B82FE */
      /* Scaled Q's: Qn_here = 2**n * Qn_above, for R(2*z) where z = hxs =
         x*x/2: */
      Q1 = -3.33333333333331316428e-02, /* BFA11111 111110F4 */
      Q2 = 1.58730158725481460165e-03,  /* 3F5A01A0 19FE5585 */
      Q3 = -7.93650757867487942473e-05, /* BF14CE19 9EAADBB7 */
      Q4 = 4.00821782732936239552e-06,  /* 3ED0CFCA 86E65239 */
      Q5 = -2.01099218183624371326e-07; /* BE8AFDB7 6E09C32D */

  static volatile double huge = 1.0e+300;

  double y, hi, lo, c, t, e, hxs, hfx, r1, twopk;
  int32_t k, xsb;
  uint32_t hx;

  GET_HIGH_WORD(hx, x);
  xsb = hx & 0x80000000; /* sign bit of x */
  hx &= 0x7FFFFFFF;      /* high word of |x| */

  /* filter out huge and non-finite argument */
  if (hx >= 0x4043687A) {   /* if |x|>=56*ln2 */
    if (hx >= 0x40862E42) { /* if |x|>=709.78... */
      if (hx >= 0x7FF00000) {
        uint32_t low;
        GET_LOW_WORD(low, x);
        if (((hx & 0xFFFFF) | low) != 0)
          return x + x; /* NaN */
        else
          return (xsb == 0) ? x : -1.0; /* exp(+-inf)={inf,-1} */
      }
      if (x > o_threshold) return huge * huge; /* overflow */
    }
    if (xsb != 0) {        /* x < -56*ln2, return -1.0 with inexact */
      if (x + tiny < 0.0)  /* raise inexact */
        return tiny - one; /* return -1 */
    }
  }

  /* argument reduction */
  if (hx > 0x3FD62E42) {   /* if  |x| > 0.5 ln2 */
    if (hx < 0x3FF0A2B2) { /* and |x| < 1.5 ln2 */
      if (xsb == 0) {
        hi = x - ln2_hi;
        lo = ln2_lo;
        k = 1;
      } else {
        hi = x + ln2_hi;
        lo = -ln2_lo;
        k = -1;
      }
    } else {
      k = invln2 * x + ((xsb == 0) ? 0.5 : -0.5);
      t = k;
      hi = x - t * ln2_hi; /* t*ln2_hi is exact here */
      lo = t * ln2_lo;
    }
    x = hi - lo;
    c = (hi - x) - lo;
  } else if (hx < 0x3C900000) { /* when |x|<2**-54, return x */
    t = huge + x;               /* return x with inexact flags when x!=0 */
    return x - (t - (huge + x));
  } else {
    k = 0;
  }

  /* x is now in primary range */
  hfx = 0.5 * x;
  hxs = x * hfx;
  r1 = one + hxs * (Q1 + hxs * (Q2 + hxs * (Q3 + hxs * (Q4 + hxs * Q5))));
  t = 3.0 - r1 * hfx;
  e = hxs * ((r1 - t) / (6.0 - x * t));
  if (k == 0) {
    return x - (x * e - hxs); /* c is 0 */
  } else {
    INSERT_WORDS(
        twopk,
        0x3FF00000 + static_cast<int32_t>(static_cast<uint32_t>(k) << 20),
        0); /* 2^k */
    e = (x * (e - c) - c);
    e -= hxs;
    if (k == -1) return 0.5 * (x - e) - 0.5;
    if (k == 1) {
      if (x < -0.25)
        return -2.0 * (e - (x + 0.5));
      else
        return one + 2.0 * (x - e);
    }
    if (k <= -2 || k > 56) { /* suffice to return exp(x)-1 */
      y = one - (e - x);
      // TODO(mvstanton): is this replacement for the hex float
      // sufficient?
      // if (k == 1024) y = y*2.0*0x1p1023;
      if (k == 1024)
        y = y * 2.0 * 8.98846567431158e+307;
      else
        y = y * twopk;
      return y - one;
    }
    t = one;
    if (k < 20) {
      SET_HIGH_WORD(t, 0x3FF00000 - (0x200000 >> k)); /* t=1-2^-k */
      y = t - (e - x);
      y = y * twopk;
    } else {
      SET_HIGH_WORD(t, ((0x3FF - k) << 20)); /* 2^-k */
      y = x - (e + t);
      y += one;
      y = y * twopk;
    }
  }
  return y;
}

double pow(double x, double y) {
  static const double
      bp[] = {1.0, 1.5},
      dp_h[] = {0.0, 5.84962487220764160156e-01},  // 0x3FE2B803, 0x40000000
      dp_l[] = {0.0, 1.35003920212974897128e-08},  // 0x3E4CFDEB, 0x43CFD006
      zero = 0.0, one = 1.0, two = 2.0,
      two53 = 9007199254740992.0,  // 0x43400000, 0x00000000
      huge = 1.0e300, tiny = 1.0e-300,
      // poly coefs for (3/2)*(log(x)-2s-2/3*s**3
      L1 = 5.99999999999994648725e-01,      // 0x3FE33333, 0x33333303
      L2 = 4.28571428578550184252e-01,      // 0x3FDB6DB6, 0xDB6FABFF
      L3 = 3.33333329818377432918e-01,      // 0x3FD55555, 0x518F264D
      L4 = 2.72728123808534006489e-01,      // 0x3FD17460, 0xA91D4101
      L5 = 2.30660745775561754067e-01,      // 0x3FCD864A, 0x93C9DB65
      L6 = 2.06975017800338417784e-01,      // 0x3FCA7E28, 0x4A454EEF
      P1 = 1.66666666666666019037e-01,      // 0x3FC55555, 0x5555553E
      P2 = -2.77777777770155933842e-03,     // 0xBF66C16C, 0x16BEBD93
      P3 = 6.61375632143793436117e-05,      // 0x3F11566A, 0xAF25DE2C
      P4 = -1.65339022054652515390e-06,     // 0xBEBBBD41, 0xC5D26BF1
      P5 = 4.13813679705723846039e-08,      // 0x3E663769, 0x72BEA4D0
      lg2 = 6.93147180559945286227e-01,     // 0x3FE62E42, 0xFEFA39EF
      lg2_h = 6.93147182464599609375e-01,   // 0x3FE62E43, 0x00000000
      lg2_l = -1.90465429995776804525e-09,  // 0xBE205C61, 0x0CA86C39
      ovt = 8.0085662595372944372e-0017,    // -(1024-log2(ovfl+.5ulp))
      cp = 9.61796693925975554329e-01,      // 0x3FEEC709, 0xDC3A03FD =2/(3ln2)
      cp_h = 9.61796700954437255859e-01,    // 0x3FEEC709, 0xE0000000 =(float)cp
      cp_l = -7.02846165095275826516e-09,   // 0xBE3E2FE0, 0x145B01F5 =tail cp_h
      ivln2 = 1.44269504088896338700e+00,   // 0x3FF71547, 0x652B82FE =1/ln2
      ivln2_h =
          1.44269502162933349609e+00,  // 0x3FF71547, 0x60000000 =24b 1/ln2
      ivln2_l =
          1.92596299112661746887e-08;  // 0x3E54AE0B, 0xF85DDF44 =1/ln2 tail

  double z, ax, z_h, z_l, p_h, p_l;
  double y1, t1, t2, r, s, t, u, v, w;
  int i, j, k, yisint, n;
  int hx, hy, ix, iy;
  unsigned lx, ly;

  EXTRACT_WORDS(hx, lx, x);
  EXTRACT_WORDS(hy, ly, y);
  ix = hx & 0x7fffffff;
  iy = hy & 0x7fffffff;

  /* y==zero: x**0 = 1 */
  if ((iy | ly) == 0) return one;

  /* +-NaN return x+y */
  if (ix > 0x7ff00000 || ((ix == 0x7ff00000) && (lx != 0)) || iy > 0x7ff00000 ||
      ((iy == 0x7ff00000) && (ly != 0))) {
    return x + y;
  }

  /* determine if y is an odd int when x < 0
   * yisint = 0 ... y is not an integer
   * yisint = 1 ... y is an odd int
   * yisint = 2 ... y is an even int
   */
  yisint = 0;
  if (hx < 0) {
    if (iy >= 0x43400000) {
      yisint = 2; /* even integer y */
    } else if (iy >= 0x3ff00000) {
      k = (iy >> 20) - 0x3ff; /* exponent */
      if (k > 20) {
        j = ly >> (52 - k);
        if ((j << (52 - k)) == static_cast<int>(ly)) yisint = 2 - (j & 1);
      } else if (ly == 0) {
        j = iy >> (20 - k);
        if ((j << (20 - k)) == iy) yisint = 2 - (j & 1);
      }
    }
  }

  /* special value of y */
  if (ly == 0) {
    if (iy == 0x7ff00000) { /* y is +-inf */
      if (((ix - 0x3ff00000) | lx) == 0) {
        return y - y;                /* inf**+-1 is NaN */
      } else if (ix >= 0x3ff00000) { /* (|x|>1)**+-inf = inf,0 */
        return (hy >= 0) ? y : zero;
      } else { /* (|x|<1)**-,+inf = inf,0 */
        return (hy < 0) ? -y : zero;
      }
    }
    if (iy == 0x3ff00000) { /* y is  +-1 */
      if (hy < 0) {
        return base::Divide(one, x);
      } else {
        return x;
      }
    }
    if (hy == 0x40000000) return x * x; /* y is  2 */
    if (hy == 0x3fe00000) {             /* y is  0.5 */
      if (hx >= 0) {                    /* x >= +0 */
        return sqrt(x);
      }
    }
  }

  ax = fabs(x);
  /* special value of x */
  if (lx == 0) {
    if (ix == 0x7ff00000 || ix == 0 || ix == 0x3ff00000) {
      z = ax;                         /*x is +-0,+-inf,+-1*/
      if (hy < 0) z = base::Divide(one, z); /* z = (1/|x|) */
      if (hx < 0) {
        if (((ix - 0x3ff00000) | yisint) == 0) {
          /* (-1)**non-int is NaN */
          z = std::numeric_limits<double>::signaling_NaN();
        } else if (yisint == 1) {
          z = -z; /* (x<0)**odd = -(|x|**odd) */
        }
      }
      return z;
    }
  }

  n = (hx >> 31) + 1;

  /* (x<0)**(non-int) is NaN */
  if ((n | yisint) == 0) {
    return std::numeric_limits<double>::signaling_NaN();
  }

  s = one; /* s (sign of result -ve**odd) = -1 else = 1 */
  if ((n | (yisint - 1)) == 0) s = -one; /* (-ve)**(odd int) */

  /* |y| is huge */
  if (iy > 0x41e00000) {   /* if |y| > 2**31 */
    if (iy > 0x43f00000) { /* if |y| > 2**64, must o/uflow */
      if (ix <= 0x3fefffff) return (hy < 0) ? huge * huge : tiny * tiny;
      if (ix >= 0x3ff00000) return (hy > 0) ? huge * huge : tiny * tiny;
    }
    /* over/underflow if x is not close to one */
    if (ix < 0x3fefffff) return (hy < 0) ? s * huge * huge : s * tiny * tiny;
    if (ix > 0x3ff00000) return (hy > 0) ? s * huge * huge : s * tiny * tiny;
    /* now |1-x| is tiny <= 2**-20, suffice to compute
       log(x) by x-x^2/2+x^3/3-x^4/4 */
    t = ax - one; /* t has 20 trailing zeros */
    w = (t * t) * (0.5 - t * (0.3333333333333333333333 - t * 0.25));
    u = ivln2_h * t; /* ivln2_h has 21 sig. bits */
    v = t * ivln2_l - w * ivln2;
    t1 = u + v;
    SET_LOW_WORD(t1, 0);
    t2 = v - (t1 - u);
  } else {
    double ss, s2, s_h, s_l, t_h, t_l;
    n = 0;
    /* take care subnormal number */
    if (ix < 0x00100000) {
      ax *= two53;
      n -= 53;
      GET_HIGH_WORD(ix, ax);
    }
    n += ((ix) >> 20) - 0x3ff;
    j = ix & 0x000fffff;
    /* determine interval */
    ix = j | 0x3ff00000; /* normalize ix */
    if (j <= 0x3988E) {
      k = 0; /* |x|<sqrt(3/2) */
    } else if (j < 0xBB67A) {
      k = 1; /* |x|<sqrt(3)   */
    } else {
      k = 0;
      n += 1;
      ix -= 0x00100000;
    }
    SET_HIGH_WORD(ax, ix);

    /* compute ss = s_h+s_l = (x-1)/(x+1) or (x-1.5)/(x+1.5) */
    u = ax - bp[k]; /* bp[0]=1.0, bp[1]=1.5 */
    v = base::Divide(one, ax + bp[k]);
    ss = u * v;
    s_h = ss;
    SET_LOW_WORD(s_h, 0);
    /* t_h=ax+bp[k] High */
    t_h = zero;
    SET_HIGH_WORD(t_h, ((ix >> 1) | 0x20000000) + 0x00080000 + (k << 18));
    t_l = ax - (t_h - bp[k]);
    s_l = v * ((u - s_h * t_h) - s_h * t_l);
    /* compute log(ax) */
    s2 = ss * ss;
    r = s2 * s2 *
        (L1 + s2 * (L2 + s2 * (L3 + s2 * (L4 + s2 * (L5 + s2 * L6)))));
    r += s_l * (s_h + ss);
    s2 = s_h * s_h;
    t_h = 3.0 + s2 + r;
    SET_LOW_WORD(t_h, 0);
    t_l = r - ((t_h - 3.0) - s2);
    /* u+v = ss*(1+...) */
    u = s_h * t_h;
    v = s_l * t_h + t_l * ss;
    /* 2/(3log2)*(ss+...) */
    p_h = u + v;
    SET_LOW_WORD(p_h, 0);
    p_l = v - (p_h - u);
    z_h = cp_h * p_h; /* cp_h+cp_l = 2/(3*log2) */
    z_l = cp_l * p_h + p_l * cp + dp_l[k];
    /* log2(ax) = (ss+..)*2/(3*log2) = n + dp_h + z_h + z_l */
    t = static_cast<double>(n);
    t1 = (((z_h + z_l) + dp_h[k]) + t);
    SET_LOW_WORD(t1, 0);
    t2 = z_l - (((t1 - t) - dp_h[k]) - z_h);
  }

  /* split up y into y1+y2 and compute (y1+y2)*(t1+t2) */
  y1 = y;
  SET_LOW_WORD(y1, 0);
  p_l = (y - y1) * t1 + y * t2;
  p_h = y1 * t1;
  z = p_l + p_h;
  EXTRACT_WORDS(j, i, z);
  if (j >= 0x40900000) {               /* z >= 1024 */
    if (((j - 0x40900000) | i) != 0) { /* if z > 1024 */
      return s * huge * huge;          /* overflow */
    } else {
      if (p_l + ovt > z - p_h) return s * huge * huge; /* overflow */
    }
  } else if ((j & 0x7fffffff) >= 0x4090cc00) { /* z <= -1075 */
    if (((j - 0xc090cc00) | i) != 0) {         /* z < -1075 */
      return s * tiny * tiny;                  /* underflow */
    } else {
      if (p_l <= z - p_h) return s * tiny * tiny; /* underflow */
    }
  }
  /*
   * compute 2**(p_h+p_l)
   */
  i = j & 0x7fffffff;
  k = (i >> 20) - 0x3ff;
  n = 0;
  if (i > 0x3fe00000) { /* if |z| > 0.5, set n = [z+0.5] */
    n = j + (0x00100000 >> (k + 1));
    k = ((n & 0x7fffffff) >> 20) - 0x3ff; /* new k for n */
    t = zero;
    SET_HIGH_WORD(t, n & ~(0x000fffff >> k));
    n = ((n & 0x000fffff) | 0x00100000) >> (20 - k);
    if (j < 0) n = -n;
    p_h -= t;
  }
  t = p_l + p_h;
  SET_LOW_WORD(t, 0);
  u = t * lg2_h;
  v = (p_l - (t - p_h)) * lg2 + t * lg2_l;
  z = u + v;
  w = v - (z - u);
  t = z * z;
  t1 = z - t * (P1 + t * (P2 + t * (P3 + t * (P4 + t * P5))));
  r = base::Divide(z * t1, (t1 - two) - (w + z * w));
  z = one - (r - z);
  GET_HIGH_WORD(j, z);
  j += static_cast<int>(static_cast<uint32_t>(n) << 20);
  if ((j >> 20) <= 0) {
    z = scalbn(z, n); /* subnormal output */
  } else {
    int tmp;
    GET_HIGH_WORD(tmp, z);
    SET_HIGH_WORD(z, tmp + static_cast<int>(static_cast<uint32_t>(n) << 20));
  }
  return s * z;
}

#undef EXTRACT_WORDS
#undef GET_HIGH_WORD
#undef GET_LOW_WORD
#undef INSERT_WORDS
#undef SET_HIGH_WORD
#undef SET_LOW_WORD
}
