generate_fcl --description DS-flash --dsconf MDC2018a --dsowner mu2e --auxinput=1:physics.filters.flashResample.fileNames:mubeam_10.txt --include JobConfig/beam/DS-flash.fcl --events-per-job 200000 --njobs 10 --run-number 10021
rm -rf DS-flash
mv 000 DS-flash
