generate_fcl --description DS-cosmic-nofield --dsconf MDC2020Dev-a --dsowner mu2e --input DS-cosmic-nofield.txt --include JobConfig/reprocess/mdc2018--mdc2020-dev/dig-DS-cosmic-nofield.fcl --merge-factor=30
for dirname in 000 001 002 003 004 005 006 007 008 009; do
 if test -d $dirname; then
  echo "found dir $dirname"
    rm -rf DS-cosmic-nofield_$dirname
    mv $dirname DS-cosmic-nofield_$dirname
  fi
done
