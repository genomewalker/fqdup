#pragma once
#include "family.hpp"
#include <memory>
#include <mutex>
#include <vector>
namespace fqdup::r8qkeep {
// A native result borrows decoded parents until all queued router/writer work
// releases it. Recycle whole scratch batches then, retaining string capacities.
class DecodedPairPool:public std::enable_shared_from_this<DecodedPairPool> {
    std::size_t rows_;
    std::mutex mutex_;
    std::vector<std::unique_ptr<std::vector<FamilyPair>>> spare_;
public:
    explicit DecodedPairPool(std::size_t rows):rows_(rows){}
    std::shared_ptr<std::vector<FamilyPair>> acquire(){
        std::unique_ptr<std::vector<FamilyPair>> value;
        {std::lock_guard<std::mutex> lock(mutex_);if(!spare_.empty()){value=std::move(spare_.back());spare_.pop_back();}}
        if(!value)value=std::make_unique<std::vector<FamilyPair>>(rows_);
        return {value.release(),[pool=shared_from_this()](std::vector<FamilyPair>* rows) noexcept {
            std::unique_ptr<std::vector<FamilyPair>> value(rows);
            try{std::lock_guard<std::mutex> lock(pool->mutex_);pool->spare_.push_back(std::move(value));}catch(...){}
        }};
    }
};
}
