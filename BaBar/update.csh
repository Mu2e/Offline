#! /bin/tcsh

foreach package ( BField BaBar BbrGeom CLHEP DetectorModel KalmanTrack \
                  ProbTools TrajGeom TrkBase difAlgebra )
  svn update ${package}
end
