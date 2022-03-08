#!/bin/bash

for fn in 151507 151508 151509 151510 151669 151670 151671 151672 151673 151674 151675 151676; do
    wget -O ../data/DLPFC/${fn}_preprocessed.h5 https://zenodo.org/record/6334774/files/${fn}_preprocessed.h5?download=1
done
