# core-typer
A Core Genome MLST (cgMLST) Typing Tool

core-typer uses the [KMA](https://bitbucket.org/genomicepidemiology/kma) aligner

## Usage

```
usage: core-typer [-h] [-v] [-t THREADS] [-p PREFIX] [--min-identity MIN_IDENTITY] [--min-coverage MIN_COVERAGE] [--R1 R1] [--R2 R2] [--scheme SCHEME]
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
  --scheme SCHEME       cgMLST scheme
  --tmpdir TMPDIR       Temporary directory (default: /tmp)
  --log-level LOG_LEVEL
                        Log level (default: info)
  --outdir OUTDIR       Output directory
```