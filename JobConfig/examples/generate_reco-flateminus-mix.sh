#!/bin/bash
generate_fcl --description reco-flateminus-mix --dsconf MDC2018h --dsowner mu2e --inputs flateminus-mix.txt --include JobConfig/reco/flateminus-mix.fcl --merge-factor=20
for dirname in 000 001 002 003 004 005 006 007 008 009; do
 if test -d $dirname; then
  echo "found dir $dirname"
    rm -rf reco-flateminus-mix_$dirname
    mv $dirname reco-flateminus-mix_$dirname
  fi
done
