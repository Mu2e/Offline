generate_fcl --description NoPrimary-mix --dsconf MDC2020DEVa --dsowner mu2e --input NP_MDC2018d.txt --include JobConfig/reprocess/mdc2018--mdc2020-dev/dig-NoPrimary-mix.fcl --merge-factor=1
for dirname in 000 001 002 003 004 005 006 007 008 009; do
 if test -d $dirname; then
  echo "found dir $dirname"
    rm -rf NoPrimary-mix_$dirname
    mv $dirname NoPrimary-mix_$dirname
  fi
done
