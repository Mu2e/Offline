universe = grid
GridResource = gt2 fnpcosg1.fnal.gov/jobmanager-condor
executable = photonBG.sh
arguments = $(Cluster) $(Process) $ENV(USER) $ENV(PWD) /mu2e/data/outstage
output = logs/photon.$(Cluster).$(Process).out
error = logs/photon.$(Cluster).$(Process).err
log = logs/photon.$(Cluster).$(Process).log
transfer_input_files = photonBG.fcl
when_to_transfer_output = ON_EXIT
notification = NEVER
queue 85
