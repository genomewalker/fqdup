# Pinned JavaScript numerical routines

Original source: https://raw.githubusercontent.com/nodejs/node/v22.20.0/deps/v8/src/base/ieee754.cc

Original SHA256: da01a54955911cfc550117988de91516c6f1aac343af28142dce619b433a67fc.
LICENSE SHA256: 81680185f0cb339021e3b73714499264de72f33b17cfdb0c14f400de6d95cf0e.

The selected exp/log/log1p/expm1/pow function bodies and word-access macros are retained in src/merge_r8qkeep/fdlibm.cpp. Only the enclosing namespace, memcpy bit-cast shim, and double division shim replace V8 build dependencies. Sun/fdlibm and V8 notices are preserved. Node/V8 are not runtime or build dependencies of the integrated command.

The routines reproduce the pinned Node22.20.0 fitter's Math operations. Compile with contraction disabled and without fast-math. Full prior byte-parity remains the acceptance gate; fixture formatting tests alone do not establish it.
