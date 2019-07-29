#!/bin/bash
generate_fcl --description reco-flatmugamma-mix --dsconf MDC2018h --dsowner mu2e --inputs flatmugamma-mix.txt --include JobConfig/reco/flatmugamma-mix.fcl --merge-factor=20
for dirname in 000 001 002 003 004 005 006 007 008 009; do
 if test -d $dirname; then
  echo "found dir $dirname"
    rm -rf reco-flatmugamma-mix_$dirname
    mv $dirname reco-flatmugamma-mix_$dirname
  fi
done
