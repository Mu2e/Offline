generate_fcl --description flatIntRMC --dsconf MDC2018a --dsowner mu2e --include JobConfig/primary/flatIntRMC.fcl --run-number 1002 --events-per-job 60000 --njobs 500
rm -rf flatIntRMC
mv 000 flatIntRMC
