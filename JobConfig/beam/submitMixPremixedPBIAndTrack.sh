#
# This file will submit a large number of jobs to mix backgrounds and track and analyze CE events
# from CD3.  PLEASE USE THIS AS AN EXAMPLE ONLY
#
mu2eart --setup=./setup.sh 
 --fcl=JobConfig/beam/MixPremixedPBIAndTrack.fcl \
 --fclinput=1:@bgHitFiles:CD3BGHitFiles.txt \
 --inputs=CD3CE.txt \
 --outstage=/mu2e/data/outstage \
 --njobs=100
