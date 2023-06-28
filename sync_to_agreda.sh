#!/usr/bin/env bash

mkdir "../AGrEDA/data/spotless/"
rsync -ihrv --include="gold_standard_*/***" --include="reference/" --include="gold_standard_*" --exclude="*" "spotless/standards/" "../AGrEDA/data/spotless/standards/"

mkdir "../AGrEDA/data/pdac/"
rsync -ihrv "pk_all.h5ad" "../AGrEDA/data/pdac/zenodo6024273/"