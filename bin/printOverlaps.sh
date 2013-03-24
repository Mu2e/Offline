#!/bin/sh
# Andrei Gaponenko, 2012

awk '/G4Exception-START/{f=1}; { if(f) print $0; }; /G4Exception-END/{f=0;}' ${1:?"Usage: $0 mu2e-log-file-with-g4.doSurfaceCheck-set-to-true"}

## for older G4 versions:
#awk '/WARNING - G4PVPlacement::CheckOverlaps/, /(overlapping|encapsulating)/ {print}' ${1:?"Usage: $0 mu2e-log-file-with-g4.doSurfaceCheck-set-to-true"}
