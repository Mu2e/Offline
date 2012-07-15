universe = grid
GridResource = gt2 fnpcosg1.fnal.gov/jobmanager-condor
executable = neutronBG.sh
arguments = $(Cluster) $(Process) $ENV(USER) $ENV(PWD) /mu2e/data/outstage
output = logs/neutron.$(Cluster).$(Process).out
error = logs/neutron.$(Cluster).$(Process).err
log = logs/neutron.$(Cluster).$(Process).log
transfer_input_files = neutronBG.fcl
when_to_transfer_output = ON_EXIT
notification = NEVER
queue 100
