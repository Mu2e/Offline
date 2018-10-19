generate_fcl --description flatproton --dsconf MDC2018a --dsowner mu2e --embed JobConfig/primary/flatproton.fcl --run-number 1002 --events-per-job 250000 --njobs 100
rm -rf flatproton
mv 000 flatproton
