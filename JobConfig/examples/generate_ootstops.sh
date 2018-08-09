generate_fcl --description ootstops --dsconf MDC2018a --dsowner mu2e --inputs ootstops.txt --include JobConfig/beam/OOTstops.fcl --merge-factor=500
rm -rf ootstops
mv 000 ootstops
