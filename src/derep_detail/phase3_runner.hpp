#pragma once
// Phase3Runner — re-parents DerepEngine::phase3_error_correct's body without
// changing a byte of it. The body lives in phase3_error_correct.inc (verbatim
// paste of the original method; guard baseline af3ac4e1). run() opens with an
// alias block that binds every DerepEngine private the body touches, so the
// included text stays textually identical while resolving member names against
// e_. Friendship (declared in DerepEngine) grants the private access.
//
// Compiler-checked completeness: a member the body uses but the alias block
// misses fails to resolve -> build error. See docs/PHASE3_REFACTOR_SPEC.md.
//
// Included INSIDE derep.cpp's anonymous namespace, AFTER the complete DerepEngine
// definition, so DerepEngine and all file-local types/functions are visible.

class Phase3Runner {
public:
    explicit Phase3Runner(DerepEngine& e) : e_(e) {}

    void run() {
        // --- alias block: 15 data members bound by reference (writes propagate)
        auto& arena_                        = e_.arena_;
        auto& qual_arena_                   = e_.qual_arena_;
        auto& errcor_                       = e_.errcor_;
        auto& errcor_absorbed_              = e_.errcor_absorbed_;
        auto& fqcl_mismatches_              = e_.fqcl_mismatches_;
        auto& fqcl_parent_chain_            = e_.fqcl_parent_chain_;
        auto& fqcl_path_                    = e_.fqcl_path_;
        auto& index_                        = e_.index_;
        auto& is_error_                     = e_.is_error_;
        auto& profile_                      = e_.profile_;
        auto& loss_bucket_overflow_drops_   = e_.loss_bucket_overflow_drops_;
        auto& loss_short_brute_evaluated_   = e_.loss_short_brute_evaluated_;
        auto& loss_short_brute_found_       = e_.loss_short_brute_found_;
        auto& loss_short_interior_skipped_  = e_.loss_short_interior_skipped_;
        auto& loss_short_too_small_skipped_ = e_.loss_short_too_small_skipped_;
        // --- write_fqcl_ is a member method; shim keeps the call site verbatim
        auto write_fqcl_ = [this](const std::vector<IndexEntry*>& seq_entry) {
            e_.write_fqcl_(seq_entry);
        };

        #include "phase3_error_correct.inc"
    }

private:
    DerepEngine& e_;
};
