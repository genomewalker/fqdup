# Performance

## Benchmarks

Tested on a 25.8 M-read ancient DNA library (`a88af16f35`, sorted FASTQ,
`~65 bp` average read length, `.gz` input):

| Mode | Unique clusters | Wall time |
|------|----------------|-----------|
| Standard | 5,582,073 | ~30 s |
| Damage-aware | 5,547,508 | ~31 s |
| Damage-aware + error-correct | 5,532,327 | ~33 s |

Hardware: 96-core server, NFS-mounted storage.

## Memory

| Component | Usage |
|-----------|-------|
| Pass 1 index | ~16 bytes × input reads |
| SeqArena (error-correct only) | ~1 byte/bp × unique reads |
| Phase 3 pair-key index | ~O(parents) |

For 400 M reads with 50% duplication rate and `--error-correct`:
- Pass 1 index: ~6.4 GB
- SeqArena: ~10–15 GB (200 M unique reads × 65 bp average)
- Total: ~17–22 GB

Without `--error-correct`, SeqArena is never allocated:
- 400 M reads: ~6.4 GB

## I/O bottleneck

Throughput is typically limited by I/O, not CPU (~40,000–50,000 reads/s on
NFS). Acceleration options:

| Flag | Speedup | Requirement |
|------|---------|-------------|
| `--isal` | 4–6× decompression | Intel ISA-L library |
| `--pigz` | 2–3× decompression | pigz in PATH |
| `--fast` (sort) | 3× sort I/O | More disk space |

On local NVMe the I/O bottleneck shrinks significantly and CPU becomes
relevant — use `--isal` in that case.

## Scaling

- `fqdup sort` scales linearly with available cores for the sort phase.
- `fqdup derep` is single-threaded (I/O bound). Phase 3 is also
  single-threaded but uses AVX2 SIMD for the Hamming checks.
- Both subcommands stream input exactly once (twice for `derep` — Pass 1
  and Pass 2) with no random access.
