generate_fcl --description flatmugamma-calo --dsconf MDC2018a --dsowner mu2e --include JobConfig/primary/flatmugamma-calo.fcl --run-number 1002 --events-per-job 500000 --njobs 10
rm -rf flatmugamma-calo
mv 000 flatmugamma-calo
