#! /bin/bash
#
# Given a directory name, descend the directory tree looking for files with
# the name *_data.root.  Write out the names of these files in a format
# that is almost suitable for pasting into a FHICL file as a list of input files.
#
# Arguments:
# 1 - the name of a directory at which the filenames are rooted.
#
# One had edit is needed: the comma following the final filename must be removed.
#
find $1 -name \*_data.root | sort |\
   awk 'BEGIN{print "inputFiles : ["} {print "    \""$1"\","} END{print "    ]"}'