generate_fcl --description flatInternalRMC --dsconf MDC2018a --dsowner mu2e --include JobConfig/primary/flatInternalRMC.fcl --run-number 1002 --events-per-job 60000 --njobs 500
rm -rf flatInternalRMC
mv 000 flatInternalRMC
