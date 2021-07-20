#! /bin/bash
#
# Copy the files from the source and build areas to create a UPS product.
#
# This script unconditionally deletes previously existing installs of the
# same product+version+qualifiers: use with care.
#

export COMPILER_CODE=e20
export ${MU2E_SETUP_BUILDOPTS}
export DEBUG_LEVEL=${build}
export PACKAGE_NAME=offline
export PACKAGE_SOURCE=${MU2E_BASE_RELEASE}
export PACKAGE_VERSION=v09_11_03

# Check that the installation directoy has been defined.
if [ "${PRODUCTS_INSTALL}" = '' ];then
    echo "The environment variable PRODUCTS_INSTALL is not set."
    echo "You must define where to install the products before sourcing this script."
    return 1
fi

# Learn if the extra products needed for the trigger are active.
# Use mu2e_artdaq_core as a proxy for the ensemble.
if [ "`ups active | grep mu2e_artdaq_core`" != "" ]; then
   haveTrigger=".trig"
else
   haveTrigger=""
fi

# There are two representations of operating system UPS flavor:
# old style, for example: Linux64bit+2.6-2.12_e7
# new style, for example: slf6.x86_64
# We need them both.
old_flavour=`ups flavor`
new_flavour=`get-directory-name subdir`

# Build the names of the directories into which we will write things
fq=${new_flavour}.${COMPILER_CODE}.${MU2E_ART_SQUALIFIER}${haveTrigger}.${DEBUG_LEVEL}
topdir=${PRODUCTS_INSTALL}/${PACKAGE_NAME}
proddir=${PRODUCTS_INSTALL}/${PACKAGE_NAME}/${PACKAGE_VERSION}
verdir=${PRODUCTS_INSTALL}/${PACKAGE_NAME}/${PACKAGE_VERSION}.version
fqdir=${PRODUCTS_INSTALL}/${PACKAGE_NAME}/${PACKAGE_VERSION}/${fq}
incdir=${PRODUCTS_INSTALL}/${PACKAGE_NAME}/${PACKAGE_VERSION}/include
cfgdir=${PRODUCTS_INSTALL}/${PACKAGE_NAME}/${PACKAGE_VERSION}/config
upsdir=${PRODUCTS_INSTALL}/${PACKAGE_NAME}/${PACKAGE_VERSION}/ups

# Make directories, if needed.
if ! [ -e ${topdir} ];then
  mkdir ${topdir}
fi

if ! [ -e ${proddir} ];then
  mkdir ${proddir}
fi

if ! [ -e ${verdir} ];then
  mkdir ${verdir}
fi

if ! [ -e ${fqdir} ];then
  mkdir ${fqdir}
fi

if ! [ -e ${cfgdir} ];then
  mkdir ${cfgdir}
fi

if ! [ -e ${incdir} ];then
  mkdir ${incdir}
fi

if ! [ -e ${upsdir} ];then
  mkdir ${upsdir}
fi

# Copy the required parts of the source directory to the installation area:

# Header files:
rsync -ar --exclude-from  ${PACKAGE_SOURCE}/ups_install/tar_exclude_for_include.txt \
    ${PACKAGE_SOURCE}  ${proddir}/include

# UPS table file
${PACKAGE_SOURCE}/ups_install/installTableFile.sh ${upsdir}/${PACKAGE_NAME}.table

# Configuration files ( .fcl, .txt and all files that will go into databases).
rsync -ar --exclude-from  ${PACKAGE_SOURCE}/ups_install/tar_exclude_for_config.txt \
    ${PACKAGE_SOURCE}  ${cfgdir}

# Libaries and binaries
rsync -ar lib ${fqdir}
rsync -ar bin ${fqdir}

# A copy of the full source
rsync -ar --exclude-from  ${PACKAGE_SOURCE}/ups_install/tar_exclude_for_source.txt \
  ${PACKAGE_SOURCE}  ${proddir}/source

# Create the ups fq files.
${PACKAGE_SOURCE}/ups_install/installFQFile.sh \
  ${verdir}/${old_flavour}_${COMPILER_CODE}_${MU2E_ART_SQUALIFIER} \
  ${COMPILER_CODE}:${MU2E_ART_SQUALIFIER}

# Repeat for the trig qualified fq files.
${PACKAGE_SOURCE}/ups_install/installFQFile.sh \
  ${verdir}/${old_flavour}_${COMPILER_CODE}_${MU2E_ART_SQUALIFIER}_trig \
  ${COMPILER_CODE}:${MU2E_ART_SQUALIFIER}:trig

unset old_flavour
unset new_flavour
unset fq
unset topdir
unset proddir
unset verdir
unset fqdir
unset incdir
unset upsdir
