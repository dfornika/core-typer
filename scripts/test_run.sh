#!/bin/bash

core-typer \
  --threads 4 \
  --scheme test_inputs/schemes/senterica/db/senterica \
  --prefix SAMN10130740 \
  --min-identity 95.0 \
  --min-coverage 90.0 \
  --R1 test_inputs/samples/SAMN10130740/SRR7903531_1.fastq.gz \
  --R2 test_inputs/samples/SAMN10130740/SRR7903531_2.fastq.gz \
  --tmpdir test_outputs/SAMN10130740/tmp \
  --outdir test_outputs/SAMN10130740 \
  --no-cleanup