#! /bin/tcsh
#
# Populate the sym link directory neeed to use the BaBar style includes.
# Execute from Offline top level directory.
#

cd BaBar
if ( -e include ) then
  cd ..
  exit
endif

mkdir include
cd include

foreach file ( BaBar BbrGeom BField CLHEP DetectorModel KalmanTrack ProbTools TrajGeom TrkBase difAlgebra )
  if !( -e $file ) then
     ln -s ../${file}/include $file
  endif
end

# Two special cases.
if !( -e ErrLogger) then
  ln -s  ../BaBar/include ErrLogger
endif
if !( -e PDT ) then
  ln -s  ../BaBar/include PDT
endif

cd ..
