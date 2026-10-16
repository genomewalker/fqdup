// Read a TDM1 damage summary with libtaph's own reader and print every field as field<TAB>value,
// in declaration order, with %.9g floats and %.17g doubles so the text round-trips exactly.
// Used by the tdm1_roundtrip test and to produce the fixture DART's reader test compares against.
#include "taph/damage_summary.hpp"

#include <cstdio>

int main(int argc, char** argv) {
    if (argc != 2) {
        std::fprintf(stderr, "usage: tdm1_dump FILE.tdm1\n");
        return 2;
    }
    const auto s = taph::read_damage_summary(argv[1]);
    if (!s) {
        std::fprintf(stderr, "tdm1_dump: cannot read TDM1 from %s\n", argv[1]);
        return 1;
    }
    std::printf("version\t%u\n", taph::DamageSummary::VERSION);
    std::printf("d_max_combined\t%.9g\n", s->d_max_combined);
    std::printf("d_max_5prime\t%.9g\n", s->d_max_5prime);
    std::printf("d_max_3prime\t%.9g\n", s->d_max_3prime);
    std::printf("lambda_5prime\t%.9g\n", s->lambda_5prime);
    std::printf("lambda_3prime\t%.9g\n", s->lambda_3prime);
    std::printf("fit_baseline_5prime\t%.9g\n", s->fit_baseline_5prime);
    std::printf("fit_baseline_3prime\t%.9g\n", s->fit_baseline_3prime);
    std::printf("terminal_shift_5prime\t%.9g\n", s->terminal_shift_5prime);
    std::printf("stop_decay_llr_5prime\t%.9g\n", s->stop_decay_llr_5prime);
    std::printf("d_max_cooc_point\t%.17g\n", s->d_max_cooc_point);
    std::printf("pi_point\t%.17g\n", s->pi_point);
    std::printf("pi_shape_lambda\t%.17g\n", s->pi_shape_lambda);
    std::printf("pi_shape_baseline\t%.17g\n", s->pi_shape_baseline);
    std::printf("pi_shape_3p_lambda\t%.17g\n", s->pi_shape_3p_lambda);
    std::printf("pi_shape_3p_baseline\t%.17g\n", s->pi_shape_3p_baseline);
    std::printf("pi_damaged_point\t%.9g\n", s->pi_damaged_point);
    std::printf("library_type\t%u\n", unsigned(s->library_type));
    std::printf("damage_validated\t%u\n", unsigned(s->damage_validated));
    std::printf("damage_artifact\t%u\n", unsigned(s->damage_artifact));
    std::printf("channel_b_valid\t%u\n", unsigned(s->channel_b_valid));
    std::printf("pi_shape_fitted\t%u\n", unsigned(s->pi_shape_fitted));
    std::printf("pi_shape_3p_fitted\t%u\n", unsigned(s->pi_shape_3p_fitted));
    std::printf("pi_state\t%u\n", unsigned(s->pi_state));
    return 0;
}
