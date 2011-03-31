#! /bin/tcsh

cd BaBar

foreach package ( BField BaBar BbrGeom CLHEP Dch DetectorModel KalmanTrack \
                  MatEnv ProbTools TrajGeom TrkBase difAlgebra )
  svn co https://opteron05.lbl.gov/svn/MU2E/${package}/trunk ${package}
end

cd ..
