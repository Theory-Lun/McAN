#!/bin/bash

../build/McAN siteMask \
  --mutation data/mut \
  --meta data/meta \
  --out out/siteMask \
  --minfreq 0.01

../build/McAN \
  --mutation data/mut \
  --meta data/meta \
  --sitemask out/siteMask \
  --outDir out \
  --oJSON \
  --oGraphML \
  --nthread 3

