#! /bin/csh -f

foreach package ( BField BaBar BbrGeom CLHEP Dch DetectorModel KalmanTrack \
                  MatEnv ProbTools TrajGeom TrkBase difAlgebra )
  svn update BaBar/${package}
end
