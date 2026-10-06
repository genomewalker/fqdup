#pragma once
#include <memory>
#include <string>
#include <vector>
namespace fqdup::r8qkeep {
struct RoutingResult { std::string routed_tsv,audit_tsv;std::vector<bool> routed; };
class RepairPipeline {
    struct Impl;
    std::shared_ptr<Impl> impl_;
public:
    RepairPipeline(const std::string& nominal_prior,const std::string& full_prior);
    RoutingResult route(const std::string& native_tsv) const;
    std::string repair(const std::string& input_tsv,int mode,double cutoff) const;
};
}
