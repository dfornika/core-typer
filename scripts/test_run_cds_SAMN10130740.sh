#!/bin/bash

core-typer \
  --threads 4 \
  --scheme test_inputs/schemes/senterica/db/senterica \
  --prefix SAMN10130740 \
  --min-identity 95.0 \
  --min-coverage 90.0 \
  --cds test_inputs/samples/SAMN10130740/GCA_004229825.1_PDT000384561.1_cds_from_genomic.fna \
  --tmpdir test_outputs/SAMN10130740/tmp_cds \
  --outdir test_outputs/SAMN10130740_cds \
  --no-cleanup
