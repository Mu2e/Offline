universe = grid
GridResource = gt2 fnpcosg1.fnal.gov/jobmanager-condor
executable = dioBG.sh
arguments = $(Cluster) $(Process) $ENV(USER) $ENV(PWD) /mu2e/data/outstage
output = logs/dio.$(Cluster).$(Process).out
error = logs/dio.$(Cluster).$(Process).err
log = logs/dio.$(Cluster).$(Process).log
transfer_input_files = dioBG.fcl
when_to_transfer_output = ON_EXIT
notification = NEVER
queue 200
