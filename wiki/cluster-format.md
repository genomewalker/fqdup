# `.fqcl` — fqdup cluster genealogy format

Persistent record of how each cluster was assembled by `fqdup derep`: parent
sequence, members, and the directed tree of single-base edits that links them.
Lets downstream tools reconstruct molecular ancestry, separate early-cycle PCR
errors from late-cycle / sequencing noise, and re-derive consensus or per-cluster
damage trajectories without re-scanning the FASTQ.

## Design goals

- **Reconstructability**: from a `.fqcl` block alone you can re-emit every
  member sequence (apply the edit chain to `parent_seq`).
- **Cluster-oriented**: one block = one cluster; random-access via a footer
  index.
- **Cheap**: ~30 B per cluster + ~6–8 B per member; orders of magnitude smaller
  than the input FASTQ.
- **Streaming-friendly**: blocks are independent and self-delimiting; the
  writer can flush without holding all clusters in memory.
- **Provenance, not consensus**: stores observed bases + edit topology;
  consensus calling is a downstream choice, not a property of the dedup.

## File layout

```
+---------------------------------------+
| MagicHeader      (16 B, fixed)        |
| MetaHeader       (variable, JSON)     |
| ClusterBlock #0  (variable)           |
| ClusterBlock #1                       |
| ...                                   |
| ClusterBlock #N                       |
| FooterIndex      (8·N + 16 B)         |
+---------------------------------------+
```

### `MagicHeader` (16 B, packed little-endian)

| offset | size | name           | value                       |
|--------|------|----------------|-----------------------------|
| 0      | 4    | magic          | `'F','Q','C','L'`           |
| 4      | 4    | version        | u32 = 1                     |
| 8      | 8    | meta_size      | u64 — bytes of MetaHeader    |

### `MetaHeader` (`meta_size` bytes, UTF-8 JSON)

```json
{
  "tool":           "fqdup",
  "tool_version":   "1.4.2",
  "input_fastq":    "merged.sorted.fq.gz",
  "n_input_reads":  113804436,
  "n_clusters":     50931408,
  "library_type":   "ds",
  "damage": {
    "d_max_5":  0.209,
    "d_max_3":  0.134,
    "lambda_5": 0.281,
    "lambda_3": 0.268,
    "mask_pos_5": 8,
    "mask_pos_3": 6
  },
  "errcor": {
    "snp_threshold":  0.20,
    "snp_min_count":  1,
    "bucket_cap":     64
  },
  "block_compression": "none"
}
```

`block_compression` ∈ `"none" | "zlib" | "zstd"` (v0.1 ships with `"none"`;
zstd is added in v0.2 once libzstd is wired into CMake).

### `ClusterBlock` (variable)

Each block is preceded by a 4-byte **length prefix** (u32, total block size in
bytes excluding itself), then the payload below. A reader can skip a block by
seeking `length` bytes forward.

| field              | type   | notes                                          |
|--------------------|--------|------------------------------------------------|
| cluster_id         | u64    | dense, 0-based, matches footer index           |
| flags              | u32    | bit0: has_member_ids, bit1: has_quality, bit2: rev_complement_used |
| n_members          | u32    | total reads collapsed into this cluster        |
| n_after_damage     | u32    | unique reads after Pass 2 damage-aware hashing |
| parent_seq_len     | u32    | bases                                          |
| parent_seq         | u8[ceil(parent_seq_len/4)] | 2-bit packed, MSB-first        |
| parent_qual_len    | u32    | 0 if `(flags & has_quality)==0`                |
| parent_qual        | u8[parent_qual_len] | raw Phred (Q0–Q60), if present    |
| n_edges            | u32    | size of the edges array                        |
| edges              | Edge[n_edges] | see below                               |
| damage_term_5      | u8[8]  | per-position C→T fraction in 5' window (×255)  |
| damage_term_3      | u8[8]  | per-position G→A fraction in 3' window (×255)  |
| n_member_ids       | u32    | 0 unless `flags & has_member_ids`              |
| member_ids         | length-prefixed string, prefix-compressed | optional        |

### `Edge` (16 B fixed)

| field          | type   | notes                                                          |
|----------------|--------|----------------------------------------------------------------|
| from_node      | u32    | parent node index in this cluster's tree (0 = parent_seq root) |
| to_node        | u32    | child node index                                               |
| pos            | u16    | position in parent_seq where the substitution occurs           |
| from_base      | u8     | 2-bit (low bits)                                               |
| to_base        | u8     | 2-bit                                                          |
| n_reads        | u32    | size of the subtree rooted at `to_node`                        |

`from_node` always precedes `to_node` in the edge order, so a reader can
reconstruct any node by walking edges sequentially and applying mutations to a
copy of `parent_seq`.

### `FooterIndex` (8·n_clusters + 16 B)

| offset (from end)         | size | name        | value                    |
|---------------------------|------|-------------|--------------------------|
| 16 + 8·n_clusters         | 8·n  | offsets[]   | u64 file offset per block |
| 16                        | 8    | n_clusters  | u64                       |
| 8                         | 4    | crc         | u32 over offsets array    |
| 4                         | 4    | magic       | `'F','Q','C','L'` (re-marker) |

A reader opens the file, seeks 16 B from the end, validates the trailing magic,
reads `n_clusters`, then seeks back `16 + 8·n_clusters` bytes to load the
offset table.

## Reconstruction algorithm

```
fn reconstruct_member(cluster, member_idx) -> bytes:
    seq = cluster.parent_seq.to_bases()       # 2-bit -> ACGT
    node = member_idx
    edits = []
    while node != 0:
        e = parent_edge_of[node]              # at most one parent per node
        edits.append((e.pos, e.to_base))
        node = e.from_node
    for (pos, base) in reverse(edits):
        seq[pos] = base
    return seq
```

`parent_edge_of` is built once per cluster by scanning `edges[]` and indexing
by `to_node` — O(n_edges) prep, O(depth) per member.

## What you can derive from `.fqcl`

| Question                                               | How                                              |
|--------------------------------------------------------|--------------------------------------------------|
| Original molecule sequence                             | `parent_seq` decode                              |
| Member set + per-member sequence                       | edge walk per node                               |
| PCR-template error vs late-PCR / sequencing            | edge `n_reads` ≫ 1 → early-cycle, propagated     |
| Per-cluster damage curve                               | `damage_term_5` / `damage_term_3`                |
| Cluster size distribution                              | iterate `n_members`                              |
| Effective unique molecule count                        | sum over clusters with `n_after_damage > 0`      |
| Variant support across clusters                        | re-emit FASTQ subset, align, count by cluster_id |

## CLI surface

```
fqdup derep ... --cluster-format clusters.fqcl
fqdup view clusters.fqcl                # cluster summary table
fqdup view clusters.fqcl --cluster 42   # ASCII tree of one cluster
fqdup view clusters.fqcl --cluster 42 --damage  # per-position damage track
fqdup export clusters.fqcl bam          # materialise members as BAM
fqdup export clusters.fqcl fastq        # one rep per cluster (or --consensus)
```

## Versioning

| version | changes                                              |
|---------|------------------------------------------------------|
| 1       | initial; raw blocks, no compression                  |
| 2 (planned) | zstd block compression + dictionary-trained `parent_seq` corpus |
| 3 (planned) | per-cluster MSA + indel-aware nodes (5-state)    |
