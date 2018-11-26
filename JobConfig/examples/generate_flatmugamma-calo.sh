generate_fcl --description flatmugamma-calo --dsconf MDC2018a --dsowner mu2e --include JobConfig/primary/flatmugamma-calo.fcl --run-number 1002 --events-per-job 250000 --njobs 200
rm -rf flatmugamma-calo
mv 000 flatmugamma-calo
