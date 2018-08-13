generate_fcl --description flashcat --dsconf MDC2018a --dsowner mu2e --inputs stashflash.txt --embed flashcat.fcl --merge-factor=200
rm -rf flashcat
mv 000 flashcat
