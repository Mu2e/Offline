#! /bin/bash
#
# Add many root files together.
#
#  $Id: addroot.sh,v 1.1 2012/07/15 16:12:11 kutschke Exp $
#  $Author: kutschke $
#  $Date: 2012/07/15 16:12:11 $
#
# Contact person Rob Kutschke
#
# Arguments:
# 1 - path to a subdidirectory
# 2 - name of the root output file
#
# Given the path (argument 1), descend from that path and find all files
# named *.root, exclude files named *_data.root. Run hadd on these
# files to form one output file. Do not add TTrees, just histograms.
#

basefile=$1
outfile=$2

inputs=`find ${basefile} -name \*.root | grep -v _data |  awk '{print $0 }'`

hadd -T $outfile $inputs
