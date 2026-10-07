#pragma once
#include <algorithm>
#include <array>
#include <cstdint>
#include <cmath>
#include <limits>
#include <string>
#include <string_view>
#include <vector>
#include <stdexcept>
namespace fqdup {
struct TerminalConstruct {
    std::string id, sequence;
    bool five_prime;
    // NaN selects the original exact policy. Learned model scores are working
    // likelihood ratios, not calibrated probabilities of technical provenance.
    double log_odds=std::numeric_limits<double>::quiet_NaN();
    bool enabled=true;
};
struct TerminalCut { size_t bases=0; std::string status="no_match", ids; };
// All coordinates refer to the original record. Matching owns no read data and
// does not allocate. The catalogue is immutable once shared with workers.
struct TerminalEvidence {
    size_t bases=0;
    const char* status="no_match";
    uint32_t first=0, last=0;
};
struct TerminalBoundaryScores {
    std::array<double,49> odds{};
    std::array<bool,49> enabled{};
    uint32_t first=0, last=0;
    bool matched=false;
};
class CompiledTerminalCatalogue {
    struct EndIndex {
        std::vector<uint32_t> offsets;
        std::vector<uint32_t> entries, short_entries;
        std::array<uint64_t,1024> present{};
    };
    std::vector<TerminalConstruct> rules_;
    std::array<EndIndex,2> ends_;
    std::vector<double> odds_;
    bool scored_=false;
    double no_cut_margin_=1e5, boundary_margin_=100;
    TerminalEvidence evaluate_scored(std::string_view s,bool five) const noexcept {
        TerminalEvidence result;
        const auto scores=boundary_scores(s,five);
        const auto& by_cut=scores.odds;
        const auto& enabled=scores.enabled;
        result.first=scores.first; result.last=scores.last;
        if(!scores.matched)return result;
        size_t best=0;double sum=1;
        for(size_t n=12;n<49;++n) {
            sum+=by_cut[n];
            if(by_cut[n]>by_cut[best]) best=n;
        }
        if(best==0) return result;
        if(by_cut[best]<no_cut_margin_) { result.status="insufficient_model_evidence";return result; }
        // Sum alternatives directly, avoiding cancellation with large odds.
        double alternatives=1;
        for(size_t n=12;n<49;++n) if(n!=best) alternatives+=by_cut[n];
        if(by_cut[best]<boundary_margin_*alternatives) { result.status="ambiguous_boundary";return result; }
        if(!enabled[best]) { result.status="unvalidated_candidate";return result; }
        result.bases=best;result.status="accepted_learned_construct";return result;
    }
    static int base(char c) noexcept {
        switch(c) { case 'A': return 0; case 'C': return 1;
                    case 'G': return 2; case 'T': return 3; default: return -1; }
    }
    // Reverse the 3-prime end, do not complement it. Thus the outer eight
    // bases are a common retrieval key for every permitted length at that end.
    static int seed(std::string_view s, bool five, bool normalize=false) noexcept {
        if(s.size()<8) return -1;
        unsigned code=0;
        for(size_t j=0;j<8;++j) {
            char c=s[five ? j : s.size()-1-j];
            if(normalize && c>='a' && c<='z') c=static_cast<char>(c-'a'+'A');
            int b=base(c);
            if(b<0) return -1;
            code=(code<<2)|static_cast<unsigned>(b);
        }
        return static_cast<int>(code);
    }
    bool matches_scored(std::string_view s,uint32_t index) const noexcept {
        const auto& r=rules_[index];if(r.sequence.size()>s.size())return false;
        size_t start=r.five_prime?0:s.size()-r.sequence.size();
        for(size_t k=0;k<r.sequence.size();++k){char c=s[start+k];if(c>='a'&&c<='z')c=static_cast<char>(c-'a'+'A');if(c!=r.sequence[k])return false;}return true;
    }
    bool matches(std::string_view s, uint32_t index) const noexcept {
        const auto& r=rules_[index];
        return r.sequence.size()<=s.size() && s.compare(
            r.five_prime ? 0 : s.size()-r.sequence.size(),r.sequence.size(),r.sequence)==0;
    }
public:
    explicit CompiledTerminalCatalogue(std::vector<TerminalConstruct> rules={},
            double no_cut_margin=1e5,double boundary_margin=100) : rules_(std::move(rules)),
            no_cut_margin_(no_cut_margin),boundary_margin_(boundary_margin) {
        if(!(no_cut_margin_>=1 && boundary_margin_>=1) || !std::isfinite(no_cut_margin_) || !std::isfinite(boundary_margin_))
            throw std::invalid_argument("invalid terminal model margins");
        for(const auto& r:rules_) if(std::isfinite(r.log_odds)) scored_=true;
        if(scored_) {
            odds_.reserve(rules_.size());
            for(const auto& r:rules_) {
                if(!std::isfinite(r.log_odds) || r.log_odds < -100 || r.log_odds > 100 || r.sequence.size()<12 || r.sequence.size()>48)
                    throw std::invalid_argument("learned configurations require finite bounded scores and 12..48 bases");
                odds_.push_back(std::exp(r.log_odds));
            }
        }
        if(rules_.size()>std::numeric_limits<uint32_t>::max())
            throw std::invalid_argument("too many terminal configurations");
        for(size_t i=0;i<rules_.size();++i) {
            const auto& r=rules_[i];
            if(r.id.empty() || r.sequence.empty() || r.sequence.size()>4096)
                throw std::invalid_argument("invalid terminal configuration length or identifier");
            for(char c:r.sequence) if(base(c)<0)
                throw std::invalid_argument("terminal configurations require A/C/G/T");
            for(size_t j=0;j<i;++j) if(r.id==rules_[j].id &&
                (r.sequence!=rules_[j].sequence || r.five_prime!=rules_[j].five_prime))
                throw std::invalid_argument("conflicting terminal configuration identifier");
        }
        if(rules_.size()<=16) return;
        for(size_t e=0;e<2;++e) {
            auto& index=ends_[e];
            index.offsets.assign(65537,0);
            for(const auto& r:rules_) if(r.five_prime==(e==0) && r.sequence.size()>=12)
                ++index.offsets[static_cast<size_t>(seed(r.sequence,r.five_prime))+1];
            for(size_t k=1;k<index.offsets.size();++k) index.offsets[k]+=index.offsets[k-1];
            index.entries.resize(index.offsets.back());
            auto cursor=index.offsets;
            for(uint32_t j=0;j<rules_.size();++j) {
                const auto& r=rules_[j];
                if(r.five_prime!=(e==0)) continue;
                if(r.sequence.size()<12) { index.short_entries.push_back(j); continue; }
                auto k=static_cast<unsigned>(seed(r.sequence,r.five_prime));
                index.entries[cursor[k]++]=j;
                index.present[k>>6]|=uint64_t{1}<<(k&63);
            }
        }
    }
    // The same alias-capped hypotheses serve the baseline and optional context
    // scorer. Disabled configurations remain competitors. No hot-path allocation.
    TerminalBoundaryScores boundary_scores(std::string_view s,bool five) const noexcept {
        TerminalBoundaryScores result;
        if(!scored_)return result;
        const auto& index=ends_[five?0:1];
        uint32_t first=0,last=static_cast<uint32_t>(rules_.size());
        if(rules_.size()>16) {
            int k=seed(s,five,true);
            if(k<0 || !(index.present[static_cast<unsigned>(k)>>6]&(uint64_t{1}<<(k&63))))return result;
            first=index.offsets[k];last=index.offsets[k+1];
            result.first=first;result.last=last;
        }
        for(uint32_t p=first;p<last;++p) {
            uint32_t j=rules_.size()<=16?p:index.entries[p];
            const auto& rule=rules_[j];
            if(rule.five_prime!=five || !matches_scored(s,j))continue;
            result.matched=true;
            size_t n=rule.sequence.size();
            result.odds[n]=std::max(result.odds[n],odds_[j]);
            result.enabled[n]=result.enabled[n]||rule.enabled;
        }
        return result;
    }
    // Offline diagnostics only. The production evaluate path stays allocation-free.
    struct BoundaryDiagnostic {
        size_t bases=0;
        double odds=0, alternatives=1;
        bool enabled=false;
        std::string ids, adjoining_context;
    };
    std::vector<BoundaryDiagnostic> explain(std::string_view s, bool five) const {
        std::array<double,49> odds{};
        std::array<bool,49> enabled{};
        std::array<std::string,49> ids{};
        for(uint32_t j=0;j<rules_.size();++j) {
            const auto& r=rules_[j];
            if(r.five_prime!=five || r.sequence.size()>48 ||
               !(scored_?matches_scored(s,j):matches(s,j))) continue;
            auto n=r.sequence.size();
            odds[n]=std::max(odds[n],scored_?odds_[j]:1.);
            enabled[n]=enabled[n]||r.enabled;
            if(!ids[n].empty()) ids[n]+=",";
            ids[n]+=r.id;
        }
        std::vector<BoundaryDiagnostic> result;
        for(size_t n=12;n<49;++n) if(!ids[n].empty()) {
            BoundaryDiagnostic d; d.bases=n;d.odds=odds[n];d.enabled=enabled[n];d.ids=ids[n];
            for(size_t k=12;k<49;++k) if(k!=n)d.alternatives+=odds[k];
            // Context is reported in original read orientation and never deleted.
            size_t count=std::min(size_t(4),s.size()-n);
            d.adjoining_context=std::string(s.substr(five?n:s.size()-n-count,count));
            result.push_back(std::move(d));
        }
        return result;
    }
    TerminalEvidence evaluate(std::string_view s, bool five) const noexcept {
        if(scored_) return evaluate_scored(s,five);
        TerminalEvidence result;
        if(rules_.size()<=16) {
            bool ambiguous=false;
            for(uint32_t j=0;j<rules_.size();++j) {
                const auto& r=rules_[j];
                if(r.five_prime!=five || !matches(s,j)) continue;
                if(r.sequence.size()<12) {
                    if(!result.bases) result.status="insufficient_fixed_match";
                    continue;
                }
                if(result.bases && result.bases!=r.sequence.size()) ambiguous=true;
                result.bases=r.sequence.size(); result.status="accepted_exact_construct";
            }
            if(ambiguous) { result.bases=0; result.status="ambiguous_boundary"; }
            return result;
        }
        const auto& index=ends_[five?0:1];
        const int k=seed(s,five,true);
        if(k>=0 && (index.present[static_cast<unsigned>(k)>>6] & (uint64_t{1}<<(k&63)))) {
            result.first=index.offsets[k]; result.last=index.offsets[k+1];
            bool ambiguous=false;
            for(uint32_t p=result.first;p<result.last;++p) {
                auto j=index.entries[p];
                if(!matches(s,j)) continue;
                size_t n=rules_[j].sequence.size();
                if(result.bases && result.bases!=n) ambiguous=true;
                result.bases=n;
                result.status="accepted_exact_construct";
            }
            if(ambiguous) { result.bases=0; result.status="ambiguous_boundary"; }
        }
        // Short explicit inputs retain the historical report-only status.
        // They do not enter the production seed table or authorize clipping.
        if(result.bases==0 && std::string_view(result.status)=="no_match")
            for(uint32_t j:index.short_entries)
                if(matches(s,j)) {
                    result.status="insufficient_fixed_match"; break;
                }
        return result;
    }
    std::string matching_ids(std::string_view s, bool five, const TerminalEvidence& evidence) const {
        std::string ids;
        if(rules_.size()<=16) {
            for(uint32_t j=0;j<rules_.size();++j)
                if(rules_[j].five_prime==five && rules_[j].sequence.size()>=12 && (scored_?matches_scored(s,j):matches(s,j)) && (!scored_ || evidence.bases==0 || rules_[j].sequence.size()==evidence.bases)) {
                    if(!ids.empty()) ids+=','; ids+=rules_[j].id;
                }
            return ids;
        }
        const auto& index=ends_[five?0:1];
        for(uint32_t p=evidence.first;p<evidence.last;++p) {
            auto j=index.entries[p];
            if(!(scored_?matches_scored(s,j):matches(s,j)) || (scored_ && evidence.bases>0 && rules_[j].sequence.size()!=evidence.bases)) continue;
            if(!ids.empty()) ids+=',';
            ids+=rules_[j].id;
        }
        return ids;
    }
    size_t bytes() const noexcept {
        size_t n=sizeof(*this);
        for(const auto& e:ends_) n+=e.offsets.capacity()*sizeof(uint32_t)+(e.entries.capacity()+e.short_entries.capacity())*sizeof(uint32_t);
        n+=rules_.capacity()*sizeof(TerminalConstruct)+odds_.capacity()*sizeof(double);
        for(const auto& r:rules_) n+=r.id.capacity()+r.sequence.capacity();
        return n;
    }
};
// Compatibility entry point for callers outside the compiled production path.
inline TerminalCut match_terminal_construct(std::string_view s,
        const std::vector<TerminalConstruct>& rules, bool five) {
    CompiledTerminalCatalogue catalogue(rules);
    auto e=catalogue.evaluate(s,five);
    return {e.bases,e.status,catalogue.matching_ids(s,five,e)};
}
inline void apply_retained_interval(std::string& seq,std::string& qual,size_t begin,size_t end) {
    if(seq.size()!=qual.size() || begin>end || end>seq.size())
        throw std::invalid_argument("invalid retained sequence/quality interval");
    if(begin==0 && end==seq.size()) return;
    seq.resize(end); qual.resize(end); seq.erase(0,begin); qual.erase(0,begin);
}
}

