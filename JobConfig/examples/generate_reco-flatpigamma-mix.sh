#!/bin/bash
generate_fcl --description reco-flatpigamma-mix --dsconf MDC2018h --dsowner mu2e --inputs flatpigamma-mix.txt --include JobConfig/reco/flatpigamma-mix.fcl --merge-factor=20
for dirname in 000 001 002 003 004 005 006 007 008 009; do
 if test -d $dirname; then
  echo "found dir $dirname"
    rm -rf reco-flatpigamma-mix_$dirname
    mv $dirname reco-flatpigamma-mix_$dirname
  fi
done
