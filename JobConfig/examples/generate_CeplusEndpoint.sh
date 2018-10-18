generate_fcl --description CeplusEndpoint --dsconf MDC2018a --dsowner mu2e --include JobConfig/primary/CeplusEndpoint.fcl --run-number 1002 --events-per-job 10000 --njobs 200
rm -rf CeplusEndpoint
mv 000 CeplusEndpoint
