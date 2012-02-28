#!/bin/sh

awk '/grid/{print > "/dev/stderr"; };{if(started) {  print $3 | "echo nz $(sort -ug|wc -l)"; print $2 | "echo ny $(sort -uug|wc -l)"; print $1 | "echo nx $(sort -uuug|wc -l)"; }}; /data/{started=1}' ${1:?"Usage: $0 bfield-text-map-file-name"}
