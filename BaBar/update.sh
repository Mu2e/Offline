#!/bin/bash
for PKG in BField BaBar BbrGeom CLHEP DetectorModel KalmanTrack MatEnv ProbTools TrajGeom TrkBase difAlgebra
do svn update BaBar/$PKG
done
