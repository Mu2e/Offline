#!/bin/bash
generate_fcl --description reco-CRY-cosmic-general-mix --dsconf MDC2018h --dsowner mu2e --inputs CRY-cosmic-general-mix.txt --include JobConfig/reco/CRY-cosmic-general-mix.fcl --merge-factor=20
for dirname in 000 001 002 003 004 005 006 007 008 009; do
 if test -d $dirname; then
  echo "found dir $dirname"
    rm -rf reco-CRY-cosmic-general-mix_$dirname
    mv $dirname reco-CRY-cosmic-general-mix_$dirname
  fi
done
