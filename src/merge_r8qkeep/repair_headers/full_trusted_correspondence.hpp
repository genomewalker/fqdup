#pragma once
#include "full_correspondence.hpp"
#include "full_complexity.hpp"
#include "self_correspondence.hpp"
#include <unordered_set>
namespace fqalign {
inline Correspondence full_trusted_correspondence(const std::string& a,const std::string& b){
 FullCorrespondence engine;const std::string mode="self";
auto r=engine.align(a,b);auto& c=r.corr;
if(c.eligible){int amin=512,amax=-1,bmin=512,bmax=-1;for(auto [i,j]:c.invariant){amin=std::min(amin,i);amax=std::max(amax,i);bmin=std::min(bmin,j);bmax=std::max(bmax,j);}
if(full_complexity::overlap_is_low_complexity(a.data()+amin,amax-amin+1)||full_complexity::overlap_is_low_complexity(b.data()+bmin,bmax-bmin+1)){c.eligible=false;c.reason="low_complexity_span";}
if(c.eligible){bool repeat=false;
if(mode=="exact20"){auto repeated=[](const std::string& str,int begin,int end){std::unordered_set<std::string> seen;for(int i=begin;i+20<=end;++i){auto k=str.substr(i,20);if(k.find('N')!=std::string::npos)continue;if(!seen.insert(k).second)return true;}return false;};repeat=repeated(a,amin,amax+1)||repeated(b,bmin,bmax+1);}
else {fqalign::SelfCorrespondence self;auto x=self.align(a,a),y=self.align(b,b);auto bad=[](const fqalign::FullResult& t){return t.min_informative>=30&&t.max_fraction<=.30&&(t.corr.paths>0);};repeat=bad(x)||bad(y);}
if(repeat){c.eligible=false;c.reason="self_repeat_ambiguity";}
}
if(c.eligible){std::vector<std::pair<int,int>> cols,used,cert;for(auto p:c.invariant)if(a[p.first]!='N'&&b[p.second]!='N')cols.push_back(p);
for(size_t begin=0;begin<cols.size();){size_t end=begin+1;while(end<cols.size()&&cols[end].first==cols[end-1].first+1&&cols[end].second==cols[end-1].second+1)++end;int total=0,before=0;for(size_t t=begin;t<end;++t)total+=a[cols[t].first]==b[cols[t].second];for(size_t t=begin;t<end;++t){auto p=cols[t];bool agree=a[p.first]==b[p.second];if(before>=5&&total-before-int(agree)>=5){used.push_back(p);if(agree)cert.push_back(p);}before+=agree;}begin=end;}
c.invariant=std::move(used);c.certified=std::move(cert);if(c.invariant.empty()){c.eligible=false;c.reason="no_anchored_columns";}}
}
return c;
}
}

