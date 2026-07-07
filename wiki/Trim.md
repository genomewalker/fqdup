# fqdup trim

Detects and removes 5'/3' adapter stubs (short adapter remnants left on collapsed/merged
reads) via hexamer frequency analysis, then clips all reads through a parallel pipeline.
User-facing subcommand (`main.cpp:68`), registered as `fqdup trim`.

Line references are to `src/trim.cpp` on `feat/merge-indel-dimer-qc` at `aa3d27f`.

## Pipeline (`trim_main` :167)

1. **Stub detection** (or manual override via `--stub5`/`--stub3`) — scan the first
   `--scan-reads` reads.
2. **Parallel clip pass** — producer → bounded clip queue → N clip workers → ordered
   output queue → writer thread (:249-337), same producer/worker/writer shape as
   `fqdup extend`.
3. **Report** — reads in/out, 5'/3' clip counts, drop count to stderr (:339-345).

## Stub detection (:210-247)

Reads up to `--scan-reads` records (default 1,000,000; 0 = all), buffering them in
`scan_buf` so the same reader is reused for the trim pass (single file open, no rescan).
While scanning:

- Reads ≥30bp (`LSD_L_MIN` :24) feed a `taph::SampleDamageProfile` via
  `taph::FrameSelector::update_sample_profile` — 5' hexamer composition and background
  damage context.
- Reads ≥12bp contribute their terminal 3' hexamer (`taph::encode_hex_at(seq, L-6)`) to a
  4096-bucket histogram (`hex3_terminal`).

`taph::detect_adapter_stubs` (libtaph, `library_interpretation_hex.cpp:130`) then compares
terminal vs. interior hexamer frequencies:

- **5' stubs**: hexamers with terminal/interior log2FC > 3.0 whose first base is not `T`
  (a `T`-leading 5' enrichment is the expected C→T deamination signal, not adapter).
  Up to 5 kept as `stubs5`.
- **3' stubs**: only computed if 5' stubs were found (gates against false positives on
  single-stranded libraries, where 3' hexamer enrichment can be genuine G→A damage).
  Same log2FC > 3.0 threshold vs. interior; hexamers ending in `A` are excluded (G→A
  deamination enriches 3'-terminal `A` on real molecules). Up to 5 kept as `stubs3`.

If neither list is populated, the run reports "No adapter stubs detected" and the clip
pass is a no-op copy.

`--stub5`/`--stub3` skip detection entirely and force a single literal 6bp stub each
(:202-209).

## Clip pass (`clip_worker` :103)

Per read, per batch (`BATCH_SZ` = 4096, :25):

- **5' clip**: if the first 6bp match any `stubs5` entry, erase them (:110-119). One clip
  max per read.
- **3' clip**: if the last 6bp match any `stubs3` entry, erase them, then repeat while a
  match remains (:120-136) — 3' stubs can chain (e.g. adapter-dimer readthrough), so this
  strips them iteratively until none match or the read drops below 12bp.
- **Length gate**: reads shorter than `--min-length` (default 15bp) after clipping have
  `seq` cleared as a drop sentinel (:137-140); the writer thread skips empty-`seq` records
  (:274).

Ordering is preserved end-to-end: batches carry a sequence id (`TrimBatch::id` :29),
workers pull from a bounded queue (`ClipQueue`, depth `2 * clip_threads`) and push results
into a `std::map`-backed reorder buffer (`OutQueue::pop_ordered` :85), so output read order
matches input regardless of which worker finished a batch first.

## CLI flags

| Flag | Effect | Default |
|---|---|---|
| `-i, --input FILE` | Input FASTQ (`.gz` or plain), required | — |
| `-o, --output FILE` | Output FASTQ (`.gz` if extension matches), required | — |
| `-p, --threads N` | Total threads; `N-1` used for clipping, 1 reserved for the writer | all cores |
| `--scan-reads N` | Reads sampled for stub detection; `0` = scan all | 1,000,000 |
| `--min-length N` | Drop trimmed reads shorter than N bp | 15 |
| `--stub5 SEQ` | Force a 5' stub, skip detection | none |
| `--stub3 SEQ` | Force a 3' stub, skip detection | none |
| `-h, --help` | Show usage | — |

## Output

Standard FASTQ (BGZF via `FastqWriter` if `-o` ends in `.gz`, plain otherwise). Records
dropped by the length gate are omitted from the output; all others are written with
their (possibly clipped) sequence and quality, in input order.
