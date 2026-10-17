#pragma once
#include "full_correction_engine.hpp"
#include "full_complexity.hpp"
#include "self_correspondence.hpp"
#include <unordered_set>
namespace fqrepair { using namespace fqalign;
inline Correspondence full_trusted_correspondence(const std::string& a,const std::string& b
#if defined(FQDUP_TRUTH_ABLATION) && FQDUP_TRUTH_ABLATION
,bool cheap=false
#endif
){
 FullCorrespondence engine;const std::string mode="self";
#if defined(FQDUP_TRUTH_ABLATION) && FQDUP_TRUTH_ABLATION
 static const std::string arm=[] {const char* p=std::getenv("FQDUP_TRUTH_ABLATION_ARM");return p?std::string(p):std::string();}();
 FullResult r;
 if(cheap&&(arm=="banded-repair"||arm=="banded-fallback")) {
   const auto g=quality_independent_geometry(a,reverse_complement(b));
   BandedCorrespondence banded;
   r.corr=banded.align(a,b,g.eligible?g.insert:-1);
   if(arm=="banded-fallback"&&!r.corr.eligible)r=engine.align(a,b);
 } else r=engine.align(a,b);
#else
 auto r=engine.align(a,b);
#endif
 auto& c=r.corr;
if(c.eligible){int amin=512,amax=-1,bmin=512,bmax=-1;for(auto [i,j]:c.invariant){amin=std::min(amin,i);amax=std::max(amax,i);bmin=std::min(bmin,j);bmax=std::max(bmax,j);}
if(full_complexity::overlap_is_low_complexity(a.data()+amin,amax-amin+1)||full_complexity::overlap_is_low_complexity(b.data()+bmin,bmax-bmin+1)){c.eligible=false;c.reason="low_complexity_span";}
if(c.eligible){bool repeat=false;
if(mode=="exact20"){auto repeated=[](const std::string& str,int begin,int end){std::unordered_set<std::string> seen;for(int i=begin;i+20<=end;++i){auto k=str.substr(i,20);if(k.find('N')!=std::string::npos)continue;if(!seen.insert(k).second)return true;}return false;};repeat=repeated(a,amin,amax+1)||repeated(b,bmin,bmax+1);}
else {fqalign::SelfCorrespondence self;auto x=self.align(a,a),y=self.align(b,b);auto bad=[](const fqalign::FullResult& t){return t.min_informative>=30&&t.max_fraction<=.30&&(t.corr.paths>0);};repeat=bad(x)||bad(y);}
if(repeat){c.eligible=false;c.reason="self_repeat_ambiguity";}
}
}
return c;
}
}

