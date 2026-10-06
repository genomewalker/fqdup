#pragma once
#include <array>
#include <istream>
#include <optional>
#include <string>
#include <vector>
namespace fqdup::r8qkeep {
constexpr int mixture_components = 66;
struct PairLikelihood {
    std::string id;
    int fold = 0;
    std::array<double, mixture_components> log_likelihood{};
};
struct MixtureFit {
    std::array<double, mixture_components> weights{};
    std::size_t training_pairs = 0;
    bool converged = false;
    int iterations = 0, polish_iterations = 0;
    double kkt_gap = 0, log_likelihood = 0;
};
// Input pairs and folds must preserve the census's first-seen ordering.
std::optional<MixtureFit> fit_mixture(const std::vector<PairLikelihood>& pairs);
std::vector<PairLikelihood> read_census_likelihood(std::istream& census);
std::string prior_text(const std::array<std::optional<MixtureFit>,2>& fits);
}
