universe = grid
GridResource = gt2 fnpcosg1.fnal.gov/jobmanager-condor
executable = protonBG.sh
arguments = $(Cluster) $(Process) $ENV(USER) $ENV(PWD) /mu2e/data/outstage
output = logs/proton.$(Cluster).$(Process).out
error = logs/proton.$(Cluster).$(Process).err
log = logs/proton.$(Cluster).$(Process).log
transfer_input_files = protonBG.fcl
when_to_transfer_output = ON_EXIT
notification = NEVER
queue 15
