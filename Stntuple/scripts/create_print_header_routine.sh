#!/usr/bin/env bash

# echo $1  
# echo .............. $2

file=$1

tag=dev_411
minor_tag=01
release=`echo $MU2ESOFT`
#-----------------------------------------------------------------------
# first check if what is being built is different from dev_243
#-----------------------------------------------------------------------
# cvs rdiff -r $tag Stntuple 2>&1 | grep -v Diffing > /tmp/$$.$USER
# ndiff=`cat /tmp/$$.$USER | wc -c`
ndiff=0

if [ -f $file ] ; then rm $file ; fi

echo  "#include <cstdio>"                                                    >  $file
echo  "#include <cstdlib>"                                                   >> $file
echo  "#include <cstring>"                                                   >> $file
echo  ""
echo  "void stntuple_print_header() {"                                       >> $file
echo  "  printf(\"stnmaker.exe built on `date` by $USER@`hostname` \n\");"   >> $file
echo  "  printf(\" using CDFSOFT release: $release \");"                     >> $file
echo  "  printf(\"and working directory $2 \n\");"                           >> $file
echo  "  printf(\"differences with ${tag}_$minor_tag: $ndiff bytes\n\");"    >> $file
echo  "}"                                                                    >> $file
echo  ""
echo  "void stntuple_get_version(char* Version, char* Text) {"               >> $file
echo  "  static char  txt[200];"                                             >> $file
echo  "  strcpy(Version,\"${tag}_${minor_tag}\");"                           >> $file
echo  "  strcpy(txt,\"stnmaker.exe ${tag}_${minor_tag} \");"                 >> $file
echo  "  strcat(txt,\"built on `date` by $USER@`hostname`\");"               >> $file
echo  "  strcat(txt,\" using CDFSOFT release: $release \");"                 >> $file
echo  "  strcat(txt,\"and working directory $2 \");"                         >> $file
echo  "  Text = txt;"                                                        >> $file
echo  "}"                                                                    >> $file


