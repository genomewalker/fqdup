#include "census.hpp"
#include "ordered_chunks.hpp"
#include "geometry_counts.hpp"
#include "base_validation.hpp"
#include "diagonal_certificate.hpp"
#include "fast_dp.hpp"
#include <string_view>
#include "string_input.hpp"
#include "fdlibm.hpp"
#include "shortest_double.hpp"
#define fqalign fqdup_r8qkeep_census_alignment
#define fqrepair fqdup_r8qkeep_census_repair
#include "native_quality.hpp"
#include "full_correspondence.hpp"
#include "full_complexity.hpp"
#include "self_correspondence.hpp"
#include "repair_correspondence.hpp"
#include <algorithm>
#include <array>
#include <iomanip>
#include <sstream>
#include <stdexcept>
namespace fqdup::r8qkeep {
namespace {
int fold(const FamilyPair& p) {std::uint64_t h=14695981039346656037ULL;for(unsigned char c:p.r[0].s+"|"+p.r[1].s){h^=c;h*=1099511628211ULL;}return h&1;}
int bin(char q) {const int edges[]={0,3,8,13,18,23,28,33,38,94},v=int(static_cast<unsigned char>(q))-33;int k=0;while(k<8&&v>=edges[k+1])++k;return k;}
fqalign::FullResult full_columns(fqalign::FullCorrespondence& engine,const std::string& a,const std::string& b) {
    auto r=engine.align(a,b);auto& c=r.corr;
    if(c.eligible) {
        int amin=512,amax=-1,bmin=512,bmax=-1;for(auto [i,j]:c.invariant){amin=std::min(amin,i);amax=std::max(amax,i);bmin=std::min(bmin,j);bmax=std::max(bmax,j);}
        if(fqalign::full_complexity::overlap_is_low_complexity(a.data()+amin,amax-amin+1)||fqalign::full_complexity::overlap_is_low_complexity(b.data()+bmin,bmax-bmin+1)){c.eligible=false;c.reason="low_complexity_span";}
        if(c.eligible){fqalign::SelfCorrespondence self;const auto x=self.align(a,a),y=self.align(b,b);const auto bad=[](const fqalign::FullResult& t){return t.min_informative>=30&&t.max_fraction<=.30&&(t.corr.paths>0);};if(bad(x)||bad(y)){c.eligible=false;c.reason="self_repeat_ambiguity";}}
        if(c.eligible) {
            std::vector<std::pair<int,int>> cols,used,cert;for(auto p:c.invariant)if(a[p.first]!='N'&&b[p.second]!='N')cols.push_back(p);
            for(std::size_t begin=0;begin<cols.size();) {
                std::size_t end=begin+1;while(end<cols.size()&&cols[end].first==cols[end-1].first+1&&cols[end].second==cols[end-1].second+1)++end;
                int total=0,before=0;for(std::size_t t=begin;t<end;++t)total+=a[cols[t].first]==b[cols[t].second];
                for(std::size_t t=begin;t<end;++t){const auto p=cols[t];const bool agree=a[p.first]==b[p.second];if(before>=5&&total-before-int(agree)>=5){used.push_back(p);if(agree)cert.push_back(p);}before+=agree;}
                begin=end;
            }
            c.invariant=std::move(used);c.certified=std::move(cert);if(c.invariant.empty()){c.eligible=false;c.reason="no_anchored_columns";}
        }
    }
    return r;
}
}
namespace {
// Both source vectors remain immutable until ordered_index_chunks drains.
// Only likelihoods/debug text are owned by each returned chunk.
struct CensusPairs {
    const std::vector<FamilyPair>& all;
    std::size_t offset,count;
    auto begin()const{return all.begin()+offset;}
    auto end()const{return begin()+count;}
    std::size_t size()const{return count;}
    bool empty()const{return count==0;}
    const FamilyPair& operator[](std::size_t i)const{return all[offset+i];}
};
CensusResult census_rows(CensusPairs pairs,NativeRecords native_rows,bool debug,int threads,std::size_t start,bool counts,WorkerPool* shared_pool=nullptr) {
    if(threads>1) {
        // Duplicate IDs share an accumulator in the frozen serial algorithm.
        // Keep that path unchanged; otherwise each pair's sums are independent.
        std::vector<std::string_view> ids;ids.reserve(pairs.size());for(const auto& p:pairs)ids.emplace_back(p.id);
        std::sort(ids.begin(),ids.end());const bool unique=std::adjacent_find(ids.begin(),ids.end())==ids.end();
        std::vector<std::string_view>().swap(ids);
        if(unique) {
            const auto count=native_rows.size();
            if(std::size_t(count)!=pairs.size())throw std::runtime_error("native census row count");
            CensusResult result;result.nominal_pairs.reserve(pairs.size());bool first=true;
            ordered_index_chunks(native_rows.size(),threads,start,[&](std::size_t offset,std::size_t n){
                const CensusPairs block{pairs.all,pairs.offset+offset-start,n};
                return census_rows(block,native_rows.slice(offset-start,n),debug,1,offset,counts);
            },[&](CensusResult part){
                const auto append=[&](std::string& target,const std::string& text){target+=first?text:text.substr(text.find('\n')+1);};
                if(debug){append(result.nominal_tsv,part.nominal_tsv);append(result.full_tsv,part.full_tsv);}
                result.proposed_edits_tsv+=part.proposed_edits_tsv;
                for(auto& p:part.nominal_pairs)result.nominal_pairs.push_back(std::move(p));
                for(auto& p:part.full_pairs)result.full_pairs.push_back(std::move(p));first=false;
            },"census",0,shared_pool);
            if(!pairs.empty())return result;
        }
    }
    LikelihoodAccumulator nominal_pairs(counts),full_pairs(counts);
    std::ostringstream nominal,full,edits;nominal.precision(17);
    nominal<<"id\tfold\ti\tk\tq1_bin\tq2_bin\te1\te2\tbase1\trc_base2\tcertified\n";
    full<<"id\tfold\ti\tk\tb1\tb2\te1\te2\ta\tb\tcert\n";
    fqalign::BandedCorrespondence banded;fqalign::FullCorrespondence full_engine;
    std::string rc;
    for(const auto& p:pairs) {
        const auto& a=p.r[0];const auto& b=p.r[1];if(a.s.size()>512||b.s.size()>512)continue;
        const auto g=fqalign::quality_independent_geometry(a.s,b.s);if(!g.eligible)continue;
        fqalign::reverse_complement_into(b.s,rc);const auto c=banded.align(a.s,rc,g.insert);if(!c.eligible)continue;
        const int f=fold(p);
        for(auto [i,j]:c.invariant){if(a.s[i]=='N'||rc[j]=='N')continue;const int k=int(b.s.size())-1-j;const auto e1=fqalign::phred_error(a.q[i]),e2=fqalign::phred_error(b.q[k]);nominal_pairs.add(p.id,f,e1,e2,a.s[i]==rc[j]);if(debug){const bool cert=std::find(c.certified.begin(),c.certified.end(),std::pair<int,int>{i,j})!=c.certified.end();nominal<<p.id<<'\t'<<f<<'\t'<<i<<'\t'<<k<<'\t'<<bin(a.q[i])<<'\t'<<bin(b.q[k])<<'\t'<<e1<<'\t'<<e2<<'\t'<<a.s[i]<<'\t'<<rc[j]<<'\t'<<cert<<'\n';}}
    }
    std::size_t row=0;
    for(const auto& v:native_rows) {
        if(row>=pairs.size())throw std::runtime_error("extra native row");const auto& p=pairs[row];const auto index=start+row++;
        if(v[0]!=p.r[0].h||v[3]!=p.r[1].h||v[1]!=p.r[0].s||v[4]!=p.r[1].s)throw std::runtime_error("native census input registry");
        if(v[6]!="unmerged")continue;
        fqalign::reverse_complement_into(v[11],rc);auto r=full_columns(full_engine,v[8],rc);const auto& c=r.corr;
        if(c.eligible){const int f=fold(p),length=int(v[11].size());for(auto [i,j]:c.invariant){const int k=length-1-j,q1=int(static_cast<unsigned char>(p.r[0].q[i]))-33,q2=int(static_cast<unsigned char>(p.r[1].q[k]))-33;const double e1=std::max(1e-10,std::min(1.0,fqdup::fdlibm::pow(10.0,-q1/10.0))),e2=std::max(1e-10,std::min(1.0,fqdup::fdlibm::pow(10.0,-q2/10.0)));const char b=fqalign::reverse_complement(std::string(1,p.r[1].s[k]))[0];full_pairs.add(p.id,f,e1,e2,p.r[0].s[i]==b);if(debug){const bool cert=std::find(c.certified.begin(),c.certified.end(),std::pair<int,int>{i,j})!=c.certified.end();full<<p.id<<'\t'<<f<<'\t'<<i<<'\t'<<k<<'\t'<<q1/5<<'\t'<<q2/5<<'\t'<<shortest_double(e1)<<'\t'<<shortest_double(e2)<<'\t'<<p.r[0].s[i]<<'\t'<<b<<'\t'<<cert<<'\n';}}}
        // Proposed edits are a debug view; repair computes its own corrections.
        // Training observations above retain the same geometry and accumulation.
        if(!debug)continue;
        const auto proposed=fqrepair::full_trusted_correspondence(v[8],rc);edits<<index<<'\t'<<proposed.eligible<<'\t'<<proposed.reason<<'\t'<<proposed.invariant.size()<<'\t';bool first=true;
        if(proposed.eligible)for(auto [i,j]:proposed.invariant){const int k=int(v[11].size())-1-j;const char a=v[8].at(i),b=rc.at(j);const int q1=int(static_cast<unsigned char>(v[9].at(i)))-33,q2=int(static_cast<unsigned char>(v[12].at(k)))-33;if(a=='N'||b=='N'||a==b||q1==q2)continue;const int mate=q1<q2?1:2,pos=mate==1?i:k,donor=mate==1?k:i;const char base=mate==1?b:fqalign::reverse_complement(std::string(1,a))[0];if(!first)edits<<';';first=false;edits<<mate<<':'<<pos<<':'<<donor<<':'<<base<<':'<<std::max(1,std::abs(q1-q2));}
        edits<<'\n';
    }
    if(row!=pairs.size())throw std::runtime_error("missing native row");
    return {debug?nominal.str():"",debug?full.str():"",edits.str(),nominal_pairs.take(),full_pairs.take()};
}
}
CensusResult census(const std::vector<FamilyPair>& pairs,const NativeResult& native,bool debug,int threads,std::size_t start,bool counts,WorkerPool* shared_pool) {
    ParsedNativeRecords parsed;
    if(native.records.empty()&&!native.native_tsv.empty())parsed=parse_native_records(native.native_tsv);
    const NativeRecords native_rows=native.records.empty()?parsed.records:native.records;
    return census_rows({pairs,0,pairs.size()},native_rows,debug,threads,start,counts,shared_pool);
}
RawConcordance raw_concordance(const std::vector<FamilyPair>& pairs,const NativeResult& native){
    if(pairs.size()!=native.records.size())throw std::runtime_error("raw concordance cardinality");
    RawConcordance result;std::ostringstream rows;rows.precision(17);
    rows<<"id\tfold\ti\tk\tq1_bin\tq2_bin\te1\te2\tbase1\trc_base2\tcertified\n";
    fqalign::BandedCorrespondence banded;std::string rc;const std::string alphabet="ACGT";
    for(std::size_t row=0;row<pairs.size();++row){
        const auto& p=pairs[row];const auto& a=p.r[0];const auto& b=p.r[1];
        if(native.records[row].outcome!="merged"||a.s.size()>512||b.s.size()>512)continue;
        const auto g=fqalign::quality_independent_geometry(a.s,b.s);
        if(!g.eligible||g.insert!=int(native.records[row].first.seq.size()))continue;
        fqalign::reverse_complement_into(b.s,rc);const auto c=banded.align(a.s,rc,g.insert);if(!c.eligible)continue;
        for(auto [i,j]:c.invariant){
            const int k=int(b.s.size())-1-j;const auto base=alphabet.find(a.s[i]);if(base==std::string::npos||alphabet.find(rc[j])==std::string::npos)continue;
            // Exactly the frozen T46 evidence probabilities, without candidate edits.
            const int q1=std::clamp(int(static_cast<unsigned char>(a.q[i]))-33,2,40),q2=std::clamp(int(static_cast<unsigned char>(b.q[k]))-33,2,40);
            const double e1=std::pow(10.,-q1/10.),e2=std::pow(10.,-q2/10.);
            const int b1=bin(char(q1+33)),b2=bin(char(q2+33));auto& cell=result.bins.cells[(b1*9+b2)*4+base];
            ++cell.n;cell.disagree+=a.s[i]!=rc[j];cell.expected+=e1+e2-4*e1*e2/3;
            const bool certified=std::find(c.certified.begin(),c.certified.end(),std::pair<int,int>{i,j})!=c.certified.end();
            rows<<p.id<<'\t'<<fold(p)<<'\t'<<i<<'\t'<<k<<'\t'<<b1<<'\t'<<b2<<'\t'<<e1<<'\t'<<e2<<'\t'<<a.s[i]<<'\t'<<rc[j]<<'\t'<<certified<<'\n';
        }
    }
    result.nominal_tsv=rows.str();return result;
}
}
