generate_fcl --description CeEndpoint --dsconf MDC2018a --dsowner mu2e --include JobConfig/primary/CeEndpoint.fcl --run-number 1002 --events-per-job 10000 --njobs 200
rm -rf CeEndpoint
mv 000 CeEndpoint
