generate_fcl --description flatInternalRPC --dsconf MDC2018a --dsowner mu2e --include JobConfig/primary/flatInternalRPC.fcl --run-number 1002 --events-per-job 60000 --njobs 500
rm -rf flatInternalRPC
mv 000 flatInternalRPC
