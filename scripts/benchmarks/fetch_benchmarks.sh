#!/bin/bash
set -euo pipefail

mkdir -p data/benchmarks
cd data/benchmarks

echo "Downloading UCSC hg38 100-way chrM..."
curl -O -s http://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz100way/maf/chrM.maf.gz

echo "Downloading UCSC cactus447way chr22 random..."
curl -O -s http://hgdownload.soe.ucsc.edu/goldenPath/hg38/cactus447way/maf/chr22_KI270735v1_random.maf.gz

echo "Downloading UCSC hg38 470-way chr22 random..."
curl -O -s http://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz470way/maf/chr22_KI270735v1_random.maf

echo "Done!"
