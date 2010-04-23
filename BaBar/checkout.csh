#! /bin/tcsh

foreach package ( BField BaBar BbrGeom CLHEP DetectorModel KalmanTrack \
                  ProbTools TrajGeom TrkBase difAlgebra )
  svn co https://opteron05.lbl.gov/svn/MU2E/${package}/trunk ${package}
end
