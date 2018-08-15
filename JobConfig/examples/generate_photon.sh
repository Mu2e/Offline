generate_fcl --description photon --dsconf MDC2018a --dsowner mu2e --embed photon.fcl --run-number 1002 --events-per-job 1000000 --njobs 10
rm -rf photon
mv 000 photon
