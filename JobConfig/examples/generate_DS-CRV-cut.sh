generate_fcl --description DS-CRV-cut --dsconf MDC2018a --dsowner mu2e --inputs DS-CRV-cat.txt --embed JobConfig/beam/DS-CRV-cut.fcl --merge-factor=1
rm -rf DS-CRV-cut
mv 000 DS-CRV-cut
