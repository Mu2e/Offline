#! /bin/bash

fqfile=$1
qualifiers=$2


cat > ${fqfile}_prof <<EOF
FILE = version
PRODUCT = ${PACKAGE_NAME}
VERSION = ${PACKAGE_VERSION}

#*************************************************
#
FLAVOR = `ups flavor`
QUALIFIERS = "${qualifiers}:prof"
  DECLARER = `whoami`
  DECLARED = `date +"%Y-%m-%d %H:%M:%S GMT" -u`
  MODIFIER = `whoami`
  MODIFIED = `date +"%Y-%m-%d %H:%M:%S GMT" -u`
  PROD_DIR = ${PACKAGE_NAME}/${PACKAGE_VERSION}
  UPS_DIR = ups
  TABLE_FILE = ${PACKAGE_NAME}.table

EOF

cat > ${fqfile}_debug <<EOF
FILE = version
PRODUCT = ${PACKAGE_NAME}
VERSION = ${PACKAGE_VERSION}

#*************************************************
#
FLAVOR = `ups flavor`
QUALIFIERS = "${qualifiers}:debug"
  DECLARER = `whoami`
  DECLARED = `date +"%Y-%m-%d %H:%M:%S GMT" -u`
  MODIFIER = `whoami`
  MODIFIED = `date +"%Y-%m-%d %H:%M:%S GMT" -u`
  PROD_DIR = ${PACKAGE_NAME}/${PACKAGE_VERSION}
  UPS_DIR = ups
  TABLE_FILE = ${PACKAGE_NAME}.table

EOF

unset fqfile
unset qualifiers
