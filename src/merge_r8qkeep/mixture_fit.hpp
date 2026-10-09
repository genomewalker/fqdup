#pragma once
#include <array>
#include <istream>
#include <optional>
#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>
namespace fqdup::r8qkeep {
constexpr int mixture_components = 66;
struct PairLikelihood {
    std::string id;
    int fold = 0;
    std::array<double, mixture_components> log_likelihood{};
    // Counts mode uses only log_likelihood[0] and these per-pair statistics.
    std::uint64_t agreements=0, disagreements=0;
    bool sufficient_statistics=false;
};
class LikelihoodAccumulator {
    std::array<std::array<double,2>,mixture_components> flat_{};
    std::unordered_map<std::string,std::size_t> slots_;
    std::vector<PairLikelihood> pairs_;
    bool counts_;
public:
    explicit LikelihoodAccumulator(bool counts=false);
    void add(const std::string& id,int fold,double e1,double e2,bool agree);
    std::vector<PairLikelihood> take(){return std::move(pairs_);}
};
struct MixtureFit {
    std::array<double, mixture_components> weights{};
    std::size_t training_pairs = 0;
    bool converged = false;
    int iterations = 0, polish_iterations = 0, newton_reentries=0;
    double kkt_gap = 0, log_likelihood = 0;
    double em_seconds=0,polish_seconds=0;
    std::uint64_t last_em_zero_cells=0,last_em_subnormal_cells=0;
    std::uint64_t exp_calls=0,exp_cache_hits=0;
    std::size_t unique_likelihoods=0;
};
// Input pairs and folds must preserve the census's first-seen ordering.
std::optional<MixtureFit> fit_mixture(const std::vector<PairLikelihood>& pairs,int threads=1);
std::optional<MixtureFit> fit_mixture_views(const std::vector<const PairLikelihood*>& pairs,int threads=1);
// One mode and training fold per call. Exact-bit cells reduce in sorted order.
// A borrowed pool must have live workers with the caller's FP environment,
// and outlive every concurrent fit call.
class WorkerPool;
std::optional<MixtureFit> fit_mixture_counts(const std::vector<const PairLikelihood*>& pairs,int threads=1,WorkerPool* shared_pool=nullptr);
std::vector<PairLikelihood> read_census_likelihood(std::istream& census,bool counts=false);
std::string prior_text(const std::array<std::optional<MixtureFit>,2>& fits);
}
