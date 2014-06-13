#!/bin/awk -f
#------------------------------------------------------------------------------
#  call: timing.awk logfile
#
#  assumed format of the log file:
#
# Module name:  No. of calls:  Mean cpu time:  Mean clk time:   Cpu StdDev:
# During Event Processing:
# ========================
# Module name:         # Calls:     Mean cpu time:        Mean clk time:   Total Cpu:
# -----------------------------------------------------------------------------------
# CT_TrackingModule	  44837	0.010204+/-0.000022   0.010342+/-0.000014   457.530
# CalibrationManager	  44837	0.000062+/-0.000004   0.000061+/-0.000000   2.800
# CalorimetryModule	  44837	0.011572+/-0.000019   0.011769+/-0.000009   518.870
# CdfEmObjectModule	  44837	0.004389+/-0.000024   0.004472+/-0.000011   196.810
# CdfMetModule		  44837	0.002857+/-0.000021   0.002882+/-0.000004   128.120
# CentralStripClusterModu 44837	0.007004+/-0.000023   0.007096+/-0.000010   314.060
# CentralStripClu-pi0reco 44837	0.003598+/-0.000023   0.003657+/-0.000008   161.310
# CesMatchingModule	  44837	0.000518+/-0.000010   0.000531+/-0.000001   23.210
# ... snip...
#  
# empty line ^ marks the end of the module list
#------------------------------------------------------------------------------

BEGIN {
  name     = "";
  nevents  = 0;
  time     = 0;
}


/During Event Processing:/ {
  ready = 1;
}

/-----------------------------------------------------------------------------------/ {
  if ( ready == 1 ) ready = 2 ;
}

{
  if ( ready == 2 ) {

    if ( $0 == "" ) {
#------------------------------------------------------------------------------
#                             empty line marks end of processing
#------------------------------------------------------------------------------
      ready = 0;
    }
    else {
#------------------------------------------------------------------------------
#                             next module
#------------------------------------------------------------------------------
      time = time + $5;
      nev  = $2;
      if (( $1 != "DHInput") && ($1 != "DHOutput") && ($1 != "FileOutput")) {
	if (nev > nevents) nevents = nev;
      }
    }
  }
}


END {
  printf ("total time, nevents, time per event: %10.3f %10i %10.5f\n",
	  time,nevents, time/nevents);
}




