generate_fcl --description DS-cosmic --dsconf MDC2018a --dsowner mu2e --include JobConfig/primary/DS-cosmic.fcl --run-number 1002 --events-per-job 20000 --njobs 50000
rm -rf DS-cosmic
mv 000 DS-cosmic
