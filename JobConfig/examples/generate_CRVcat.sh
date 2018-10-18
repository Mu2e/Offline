generate_fcl  --description PS-CRVcat --dsconf MDC2018a --dsowner mu2e --embed PS-CRVcat.fcl --run-number 1002 --auxinput=100:physics.filters.crvFilter.fileNames:PS-CRVfiles.txt --events-per-job 19007 --njobs 25
rm -rf PS-CRVcat
mv 000 PS-CRVcat
generate_fcl  --description TS-CRVcat --dsconf MDC2018a --dsowner mu2e --embed TS-CRVcat.fcl --run-number 1002 --auxinput=200:physics.filters.crvFilter.fileNames:TS-CRVfiles.txt --events-per-job 487 --njobs 25
rm -rf TS-CRVcat
mv 000 TS-CRVcat
generate_fcl  --description DS-CRVcat --dsconf MDC2018a --dsowner mu2e --embed DS-CRVcat.fcl --run-number 1002 --auxinput=200:physics.filters.crvFilter.fileNames:DS-CRVfiles.txt --events-per-job 56251 --njobs 25
rm -rf DS-CRVcat
mv 000 DS-CRVcat
