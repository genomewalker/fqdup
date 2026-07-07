# fqdup sort

External merge sort of a FASTQ file by read ID: parallel chunk creation + sort,
then a k-way streaming merge. Required before `derep_pairs` / `derep`, which
assume sorted input (`main.cpp:43-48`).

Line references are to `src/sort.cpp` on `feat/merge-indel-dimer-qc` at `aa3d27f`.

## Pipeline

**Phase 1 — chunk creation** (`create_chunks` :830). A single `SortReaderBridge`
(:58, dispatches to rapidgzip/ISA-L/zlib via `make_fastq_reader`) reads records
into a pool of arena-backed chunk buffers. `num_writers = threads/2` (min 2)
worker threads pull filled buffers off a queue, sort each in place, and write
it to an uncompressed temp file (`write_sorted_chunk` :1059). The reader blocks
on backpressure until a buffer returns to the pool.

**Phase 2 — k-way merge** (`merge_chunks` :1133). Opens every sorted chunk file
via `make_fastq_reader`, merges with a min-heap keyed on read ID, and writes the
output (compressed if `-o` ends in `.gz`). Ties break by chunk index so equal
keys from different chunks always emerge in chunk order — deterministic output
regardless of which writer thread finished a chunk first (:1176-1191).

## Sort key

`trim_id` (:209): strips a leading `@` and everything from the first
space/tab.

- Default — lexicographic byte comparison of the trimmed ID.
- `-N` (natural order) — key is `(alpha prefix, numeric suffix)`; the trailing
  digit run is parsed into a `uint64_t` so `read_2` sorts before `read_10`
  (`write_sorted_chunk` :1077-1100, `merge_chunks` `make_nat_key` :1159-1167).

Both chunk sort and merge use `stable_sort` / index tie-break, so equal-key
reads keep their original relative order end to end.

## Chunk sizing

`--max-memory` drives `chunk_size_bytes_` via `adjust_chunk_size_for_input`
(:778): estimates uncompressed size assuming a 3.5x gzip ratio, then targets
`threads/2` chunks below 10GB estimated size, `threads` between 10-100GB, and
`threads*2` above 100GB, clamped to `[2GB, max_memory/20]` per chunk.
`create_chunks` additionally caps `chunk_size_bytes_` to
`max_memory / (num_writers+2) / 3` so the preallocated buffer pool stays inside
the memory budget (:843-852).

## Compression backends

Read: ISA-L igzip if built `HAVE_ISAL`, else zlib — same for write. `sort.cpp`
cannot link rapidgzip directly (`fastq_io_backend.cpp` owns those symbols;
ODR constraint, :741-742); rapidgzip reads only happen indirectly through
`SortReaderBridge` → `make_fastq_reader` (:833-834).

## CLI flags

(`print_usage` :1248-1267)

| Flag | Effect | Default |
|---|---|---|
| `-i FILE` | Input FASTQ (`.gz` supported) | required |
| `-o FILE` | Output FASTQ (`.gz` compresses) | required |
| `--max-memory SIZE` | Memory budget, e.g. `64G` | required |
| `-p THREADS` | Sorting threads | `hardware_concurrency()` |
| `-t DIR` | Temp directory | `.` |
| `-N` | Natural-order sort (numeric suffix) | off (lexicographic) |
| `--fast` | Uncompressed intermediate chunks | off |

## User-facing

Yes — top-level subcommand, first step of the documented workflow
(`main.cpp:43-48`).
