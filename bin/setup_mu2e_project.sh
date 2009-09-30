#! /bin/sh
#
# $Id: setup_mu2e_project.sh,v 1.1 2009/09/30 22:57:47 kutschke Exp $
# $Author: kutschke $
# $Date: 2009/09/30 22:57:47 $
#
# Original author Rob Kutschke
#
# Initialize the current directory tree as a Mu2e project.
#  - add the local lib area to the LD_LIBRARY_PATH
# 

if [ "`basename $0 2>/dev/null`" = "setup_mu2e_user.sh" ];then
    echo "you should be sourcing this file"; exit
fi
if [ "${MU2E_HOME-}" = '' ];then
    echo "MU2E_HOME is not set; bin/setup_mu2e.sh (from a base release)"
    echo "needs to be sourced"
    exit
fi
bin_dir=`dirname ${BASH_SOURCE}`   # assume file is in bin subdir
bin_dir=`cd $bin_dir >/dev/null 2>&1 && echo $PWD`
user_root=`dirname $bin_dir`
add_to_var $user_root/lib  LD_LIBRARY_PATH

unset bin_dir user_root
