#! /bin/sh
#
# $Id: setup_mu2e_project.sh,v 1.9 2012/09/24 19:00:42 gandr Exp $
# $Author: gandr $
# $Date: 2012/09/24 19:00:42 $
#
# Original author Rob Kutschke
#
# Initialize the current directory tree as a the root of a framework based project.
#  - add the local lib subdirectory to the LD_LIBRARY_PATH
#  - add the local Config subdirectory to the PYTHONPATH

function dropit_from_var() {
  local path="${1}"; shift
  local var
  if [[ "$1" != "-"* ]]; then var="${1}"; shift; fi
  var=${var:-PATH}
  local new_path="$(eval dropit "${@}" -e -p\"$`echo $var`\" \"\$path\")"
  eval export "$var"=\""${new_path}"\"
}

function add_to_var() {
  dropit_from_var "${@}" -s -f
}

function append_to_var() {
  dropit_from_var "${@}" -s
}

function rm_from_var() {
  dropit_from_var "${@}"
}

if [ "`basename $0 2>/dev/null`" = "setup_mu2e_project.sh" ];then
    echo "You should be sourcing this file"; exit
fi

if [ "${MU2E_BASE_RELEASE}" = '' ];then
    echo "MU2E_BASE_RELEASE is not set; "
    echo "You need to setup a base release of Mu2e Offline sourcing this file."
    return 21
fi

# Assume that this file lives one level below the root of he base release.
# So step up one level of path to define the base release.
bin_dir=`dirname ${BASH_SOURCE}`
bin_dir=`cd $bin_dir >/dev/null 2>&1 && /bin/pwd`
user_root=`dirname $bin_dir`

add_to_var $user_root/lib    LD_LIBRARY_PATH
add_to_var $user_root/scripts/build/python PYTHONPATH
add_to_var $user_root/bin    PATH

unset bin_dir user_root
