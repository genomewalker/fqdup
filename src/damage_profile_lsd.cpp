// LSD helpers: compute_lsd_edges, make_lsd_classify_params, lsd_accumulate*.
//
// Moved to libtaph (taph/lsd_accumulator.hpp, taph::* namespace) so DART and
// any other libtaph consumer share the identical per-read accumulation and
// LLR-classification math instead of reimplementing it. fqdup now only
// aliases the taph:: names (see include/fqdup/damage_profile.hpp) and keeps
// a thin make_lsd_classify_params(const DamageProfile&) wrapper (inline, in
// the header) that adapts fqdup's own bulk-fit type to the tool-independent
// taph::BulkDamageForClassify view. Nothing left to define here.

#include "fqdup/damage_profile.hpp"
