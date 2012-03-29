#!/bin/sh
# Andrei Gaponenko, 2012

awk '/WARNING - G4PVPlacement::CheckOverlaps/, /(overlapping|encapsulating)/ {print}' ${1:?"Usage: $0 mu2e-log-file-with-g4.doSurfaceCheck-set-to-true"}
