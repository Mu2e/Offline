#!/bin/awk -f
#
# Extract the primary trajectory from log file of a
#
#    mu2e -c ExtinctionMonitorFNAL/test/trajectoryIn.fcl > out.log
#
# job.  Usage:
#
#    ./ExtinctionMonitorFNAL/test/extractTrajectory.awk  out.log | tac
#
# A. Gaponenko, 2017

/^Pre:/ {
    printf( $3"\t"$4"\t"$5"\t%8.6f\t%8.6f\n", -$6/$8, -$7/$8);
}
