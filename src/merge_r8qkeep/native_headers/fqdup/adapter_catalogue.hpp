#pragma once
// Sequence constants shared with merge. Geometry and acceptance remain caller policy.
namespace fqdup::adapters {
inline constexpr const char* catalogue_version = "fqdup-terminal-v1";
inline constexpr const char* truseq_3p = "AGATCGGAAGAG";
inline constexpr const char* truseq_5p = "CTCTTCCGATCT";
inline constexpr const char* universal = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
inline constexpr const char* p7_index = "ATCTCGTATGCCGTCTTCTGCTTG";
inline constexpr const char* p5_read2 = "GTGTAGATCTCGGTGGTCGCCGTATCATT";
inline constexpr const char* ss_readthrough = "GGAAGAGCGTCGTGTAGGGAAAGAGTGT";
inline constexpr const char* ds_readthrough = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";
}

