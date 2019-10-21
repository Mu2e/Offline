# You must have previously setup mu2efiletools for this script to run
dsconf=MDC2018a
for mixin in \
oot-CRV-cat-cat \
neutron-CRV-cat \
dio-CRV-cat \
photon-CRV-cat \
PS-CRV-cut-cat \
TS-CRV-cat-cat \
DS-CRV-cut-cat \
oot-TrkCal-cat-cat \
neutron-TrkCal-cat \
dio-TrkCal-cat \
photon-TrkCal-cat \
DS-flash-TrkCal-cut-cat \
proton-TrkCal \
deuteron-TrkCal ; do
  outfile="$mixin.txt"
  dataset="sim.mu2e.$mixin.$dsconf.art"
  echo "Writing files for $dataset to $outfile"
  mu2eDatasetFileList $dataset > $outfile
done
