#!/bin/bash

core-typer \
  --threads 4 \
  --scheme test_inputs/schemes/senterica/db/senterica \
  --prefix SAMN11603843 \
  --min-identity 95.0 \
  --min-coverage 90.0 \
  --R1 test_inputs/samples/SAMN11603843/SRR9028498_1.fastq.gz \
  --R2 test_inputs/samples/SAMN11603843/SRR9028498_1.fastq.gz \
  --tmpdir test_outputs/SAMN11603843/tmp_reads \
  --outdir test_outputs/SAMN11603843_reads \
  --no-cleanup
