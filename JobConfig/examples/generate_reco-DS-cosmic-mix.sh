#!/bin/bash
generate_fcl --description reco-DS-cosmic-mix --dsconf MDC2018i --dsowner mu2e --inputs DS-cosmic-mix-cat.txt --include JobConfig/reco/DS-cosmic-mix.fcl --merge-factor=6
for dirname in 000 001 002 003 004 005 006 007 008 009; do
 if test -d $dirname; then
  echo "found dir $dirname"
    rm -rf reco-DS-cosmic-mix_$dirname
    mv $dirname reco-DS-cosmic-mix_$dirname
  fi
done
