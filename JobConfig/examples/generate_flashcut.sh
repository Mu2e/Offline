generate_fcl --description DS-flash-TrkCal-cut --dsconf MDC2018a --dsowner mu2e --inputs flash-cat.txt --embed JobConfig/beam/DS-flash-TrkCal-cut.fcl --merge-factor=1
rm -rf DS-flash-TrkCal-cut
mv 000 DS-flash-TrkCal-cut
