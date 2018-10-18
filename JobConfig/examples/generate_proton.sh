generate_fcl --description proton --dsconf MDC2018a --dsowner mu2e --include JobConfig/primary/proton.fcl --run-number 1002 --events-per-job 1000000 --njobs 10
rm -rf proton
mv 000 proton
