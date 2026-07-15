# core-typer

[![tests](https://github.com/dfornika/core-typer/actions/workflows/tests.yml/badge.svg)](https://github.com/dfornika/core-typer/actions/workflows/tests.yml)

A Core Genome MLST (cgMLST) Typing Tool

core-typer types samples from paired-end reads via [KMA](https://bitbucket.org/genomicepidemiology/kma),
or from an assembly/contigs FASTA via `blastn` (NCBI BLAST+); its allele-calling
approach is heavily influenced by CGE's
[cgMLSTFinder](https://bitbucket.org/genomicepidemiology/cgmlstfinder).

## Requirements

- Python >= 3.9
- For read input (`--R1`/`--R2`): [kma](https://bitbucket.org/genomicepidemiology/kma) on `PATH` (not a Python/pip dependency - install separately, e.g. `conda install -c bioconda kma`), and a cgMLST scheme indexed with `kma_index` (produces `<scheme>.name`, `<scheme>.comp.b`, etc.)
- For assembly input (`--assembly`): NCBI `blastn`/`makeblastdb` on `PATH` (e.g. `conda install -c bioconda blast`), and a scheme indexed with `makeblastdb -dbtype nucl` (produces `<scheme>.nin`, `<scheme>.nsq`, etc.) - see "Assembly input" below for an important naming requirement
- The two indexes can coexist under the same `--scheme` prefix if you want to support both input modes against one scheme

## Installation

```
pip install -e .
```

## Usage

```
usage: core-typer [-h] [-v] [-t THREADS] [-p PREFIX] [--min-identity MIN_IDENTITY] [--min-coverage MIN_COVERAGE]
                  [--R1 R1] [--R2 R2] [--assembly ASSEMBLY] [--scheme SCHEME]
                  [--tmpdir TMPDIR] [--log-level LOG_LEVEL] [--outdir OUTDIR]

A cgMLST Typing Tool

options:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -t THREADS, --threads THREADS
                        Number of CPU threads to use (default: 1)
  -p PREFIX, --prefix PREFIX
                        Prefix for output files (default: taken from R1 fastq file name)
  --min-identity MIN_IDENTITY
                        Minimum percent identity (default: 100.0)
  --min-coverage MIN_COVERAGE
                        Minimum percent coverage (default: 100.0)
  --R1 R1               Read 1
  --R2 R2               Read 2
  --assembly ASSEMBLY   Assembly/contigs fasta (alternative to --R1/--R2)
  --scheme SCHEME       cgMLST scheme
  --tmpdir TMPDIR       Temporary directory (default: ./tmp)
  --no-cleanup          Do not cleanup temporary directory
  --log-level LOG_LEVEL
                        Log level (default: info)
  --outdir OUTDIR       Output directory
```

Provide either `--R1`/`--R2` (reads, via `kma`) or `--assembly` (contigs, via `blastn`) - not both.

### Assembly input

`--assembly` requires the scheme's combined FASTA to be kept alongside the
blast db, named `<scheme>.fasta` (i.e. use the same path for both `-in` and
as the `-out` prefix's basename):

```
makeblastdb -in scheme_db.fasta -dbtype nucl -out scheme_db
core-typer --assembly assembly.fasta --scheme scheme_db --outdir out
```

This requirement exists because `blastdbcmd -entry all` (the obvious way to
list every template already in a blast db) turned out to be unreliable
across blast+ versions for non-accession-style headers like ours - reading
the original FASTA directly is simpler and doesn't depend on it.

kma's own assembly/contigs preset (`-asm`) was tried first and rejected:
empirically, against a synthetic scheme it missed 10-15% of loci that were
exact matches to a catalogued allele, for reasons that weren't resolved by
relaxing kma's `-p` or `-mrs`. `blastn` found every exact match reliably in
the same tests, so assembly input goes through `blastn` instead of `kma`.

Two known differences from `--R1`/`--R2` mode:
- `depth` is always blank - there's no read-depth concept for an already-assembled contig.
- `novel_allele_hash`/`novel_alleles.fsa` extraction isn't implemented yet for assembly input (it relies on kma's `-ef` consensus output, which has no `blastn` equivalent) - these are always blank/empty for `--assembly` runs.

`possible_multicopy_loci.csv`/`num_hits` still work for assembly input:
`blastn`, unlike kma's competitive read mapping, reports every hit above its
significance threshold, including several similar catalogued alleles all
matching the same genomic region - core-typer collapses hits that overlap
substantially on the query down to the single best-scoring one before
calling, so a locus with several similar catalogued alleles isn't mistaken
for a multi-copy locus.

## Output files

Calling an allele requires `query_identity` (mismatches within the aligned
region) and `template_coverage` (how much of the template was assembled) to
both meet `--min-identity`/`--min-coverage` - not `template_identity`, which
conflates the two into one number and can't tell you which one failed. This
matches CGE cgMLSTFinder's approach.

- `allele_calls.csv` - one row per locus in the scheme:
  - `call_status`: `called` (best hit met `--min-identity`/`--min-coverage`), `divergent` (a hit was found but didn't meet the thresholds), or `no_hit` (no reads mapped to this locus at all)
  - `allele_id`: the called allele id, or `-` if not `called`
  - `closest_allele_id` / `closest_allele_hash`: the allele id (and, for hash-prepared schemes, md5 hash - see below) of the best-scoring hit regardless of call_status - for `divergent` loci this is the nearest known allele; blank for `no_hit`
  - `novel_allele_hash`: for a `divergent` locus with full template coverage and a clean (unambiguous) consensus, the md5 hash of the sample's own observed sequence at that locus - a candidate for submission as a genuinely new allele. Also written out as FASTA in `novel_alleles.fsa`. Blank otherwise.
  - `query_identity`, `template_identity`, `template_coverage`, `depth`, `score`, `template_length`: from the best-scoring hit; blank for `no_hit`
  - `num_hits`: total number of candidate templates found for this locus (usually 1). Only the single best-scoring hit is ever called - a value `> 1` means more than one candidate was found and the others were discarded; see `possible_multicopy_loci.csv`.
- `novel_alleles.fsa` - FASTA of the observed consensus sequence for each `divergent`+full-coverage+clean locus, headers carrying the locus id and its `novel_allele_hash`. Empty if there are none. Reads input only (see "Assembly input" below).
- `possible_multicopy_loci.csv` - one row per candidate hit, for any locus with more than one (`num_hits > 1` in `allele_calls.csv`). A genuine second gene copy can result in more than one hit being found for a locus instead of it being silently merged into the first - this is evidence of a *possible* duplication/multi-copy locus, not confirmation (core-typer doesn't yet combine or otherwise specially handle multi-copy loci - it still just calls the single best hit in `allele_calls.csv`). Empty if there's no such evidence.
- `allele_profile.csv` - a two-row locus/allele_id matrix (header + one sample row) for downstream cgMLST comparison tools; always `-` for anything not `called`, regardless of `divergent` vs `no_hit`
- `qc.csv` - `num_loci`, `num_called_alleles`, `num_divergent`, `num_no_hit`, `num_possible_multicopy_loci`, `percent_called`, `mean_depth`, `stdev_depth` (`mean_depth`/`stdev_depth` are computed over loci with any hit - `called` or `divergent` - since `no_hit` loci have no depth measurement)

### Hash-based allele identifiers

`scripts/prepare_cgmlst_scheme.py` rewrites a cgmlst.org-style per-locus
FASTA download into one combined FASTA (suitable as input for both
`kma_index` and `makeblastdb`), with headers of the form
`<locus_id>_<allele_id>_<md5_hash>` (the md5 is of the allele's own
sequence). This lets `closest_allele_hash` show up in `allele_calls.csv` for
called/divergent loci, so an allele can be cross-referenced by its content
hash regardless of ID/numbering differences between schemes or tools.
core-typer's parsers detect this hash suffix automatically and fall back
cleanly to plain `<locus_id>_<allele_id>` schemes that don't have one.

## Testing

GitHub Actions (`.github/workflows/tests.yml`) runs the full suite - including
the `kma`/`wgsim` integration tests, installed via bioconda in CI - on every
push/PR to `main`.

Unit tests (no `kma`/`blastn` required) cover the parsers, allele-calling, and QC logic:

```
pytest tests/ -v
```

`tests/test_integration.py` (reads via `kma`) and
`tests/test_integration_assembly.py` (assembly via `blastn`) run core-typer
end-to-end against a small synthetic cgMLST scheme with a known
ground-truth allele profile. Each is skipped automatically unless its
required tools (`kma`/`kma_index`/`wgsim`, or `blastn`/`makeblastdb`,
respectively) are on `PATH`.

### Validation / threshold exploration

`scripts/synthetic_scheme.py` generates a tiny synthetic scheme (random loci,
a few SNP-derived alleles each) plus a "true genome" FASTA with a known
allele profile - no real scheme data needed. `scripts/simulate_reads.py`
wraps `wgsim` to simulate paired-end reads from that genome at a chosen depth
and error rate. Together they're useful for characterizing how
`--min-identity`/`--min-coverage` and sequencing depth affect `percent_called`
and failure modes, independent of any real scheme:

```
python3 scripts/synthetic_scheme.py --outdir /tmp/scheme --num-loci 20 --alleles-per-locus 3
kma_index -i /tmp/scheme/scheme.fasta -o /tmp/scheme/scheme_db
python3 scripts/simulate_reads.py --genome-fasta /tmp/scheme/true_genome.fasta --outdir /tmp/reads --depth 30 --error-rate 0.02
core-typer --R1 /tmp/reads/R1.fastq --R2 /tmp/reads/R2.fastq --scheme /tmp/scheme/scheme_db --outdir /tmp/out
```