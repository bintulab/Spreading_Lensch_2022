#!/bin/bash
for f in *_dedup.scaled.bedgraph
do
 echo "Processing $f"
 grep 'no_arms_puro-pA-9xTO-pEF-H2B-cit-5kbspacer-H2B-mCherry' $f > ${f%_edited.dedup.scaled.bedgraph}_reporteronly_dedup.scaled.bedgraph
done
