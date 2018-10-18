generate_fcl --description tgtstops --dsconf MDC2018a --dsowner mu2e --inputs tstops.txt --include JobConfig/beam/TGTstops.fcl --merge-factor=500
rm -rf tgtstops
mv 000 tgtstops
