generate_fcl --description PS-CRV-cut --dsconf MDC2018a --dsowner mu2e --inputs PS-CRV-cat.txt --embed JobConfig/beam/PS-CRV-cut.fcl --merge-factor=1
rm -rf PS-CRV-cut
mv 000 PS-CRV-cut
