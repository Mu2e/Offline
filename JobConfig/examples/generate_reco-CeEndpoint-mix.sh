#!/bin/bash
generate_fcl --description reco-CeEndpoint-mix --dsconf MDC2018h --dsowner mu2e --inputs CeEndpoint-mix.txt --include JobConfig/reco/CeEndpoint-mix.fcl --merge-factor=20
for dirname in 000 001 002 003 004 005 006 007 008 009; do
 if test -d $dirname; then
  echo "found dir $dirname"
    rm -rf reco-CeEndpoint-mix_$dirname
    mv $dirname reco-CeEndpoint-mix_$dirname
  fi
done
