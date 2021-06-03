#! /bin/bash
#
# Create the table file in the products directory
#
# Arguments:
#  1 - path to the table file.
#
# This script overwrites previously existing installs of the table file.
#

# Arguments
destination_file=$1

# Discover the versions of the products on which BTrk depends
art_ver=`ups active | awk '$1 == "art" {print $2}'`
heppdt_ver=`ups active | awk '$1 == "heppdt" {print $2}'`
btrk_ver=`ups active | awk '$1 == "BTrk" {print $2}'`
kinkal_ver=`ups active | awk '$1 == "KinKal" {print $2}'`
xerces_ver=`ups active | awk '$1 == "xerces_c" {print $2}'`
mu2e_artdaq_core_ver=`ups active | awk '$1 == "mu2e_artdaq_core" {print $2}'`

# A table file needs the qualifiers string in three formats:
#   - colon delimited
#   - dot deliminted
#   - colon delimited with plus signs.
# Define them all here and use them in the right spots in the table file.

# Value of the "Qualifiers " keyword, excluding trig and debug/prof/opt
qualifiers_value=${COMPILER_CODE}:${MU2E_ART_SQUALIFIER}

# Value that will be put into OFFLINE_FQ, excluding trig and debug/prof/opt
offline_fq_value=${COMPILER_CODE}.${MU2E_ART_SQUALIFIER}

# Value that will be put into MU2E_UPS_QUALFIERS, which is used
# to specify qualifiers for many of the ups products on which
# offline depends.
# Format: qual1:+qual2:+qual3
#         No leading + or colon on the first qualifer
mu2e_ups_qualifiers=${COMPILER_CODE}

# Write the table file in place
cat > ${destination_file} <<EOF
File    = table
Product = ${PACKAGE_NAME}

#*************************************************
# Starting Group definition
Group:

Flavor     = ANY
Qualifiers = "${qualifiers_value}:debug"

  Action = GetFQDir
      envSet( \${UPS_PROD_NAME_UC}_FS, "" )
      setupRequired( cetpkgsupport )
      execute( "get-directory-name subdir", NO_UPS_ENV, \${UPS_PROD_NAME_UC}_FS )
      envSet (OFFLINE_FQ, \${\${UPS_PROD_NAME_UC}_FS}.${offline_fq_value}.debug)
      envSet (MU2E_UPS_QUALIFIERS, +${mu2e_ups_qualifiers}:+debug )
      setupRequired( art  ${art_ver} -q +${mu2e_ups_qualifiers}:+debug )
      setupRequired( BTrk  ${btrk_ver} -q +${mu2e_ups_qualifiers}:+debug:+p383b )
      setupRequired( KinKal  ${kinkal_ver} -q +${mu2e_ups_qualifiers}:+debug:+p383b )
      setupRequired( heppdt  ${heppdt_ver} -q +${mu2e_ups_qualifiers}:+debug )
      setupRequired( xerces_c  ${xerces_ver} -q +${mu2e_ups_qualifiers}:+debug )

Flavor     = ANY
Qualifiers = "${qualifiers_value}:prof"

  Action = GetFQDir
      envSet( \${UPS_PROD_NAME_UC}_FS, "" )
      setupRequired( cetpkgsupport )
      execute( "get-directory-name subdir", NO_UPS_ENV, \${UPS_PROD_NAME_UC}_FS )
      envSet (OFFLINE_FQ, \${\${UPS_PROD_NAME_UC}_FS}.${offline_fq_value}.prof)
      envSet (MU2E_UPS_QUALIFIERS, +\${mu2e_ups_qualifiers}:+prof )
      setupRequired( art  ${art_ver} -q +${mu2e_ups_qualifiers}:+prof )
      setupRequired( BTrk  ${btrk_ver} -q +${mu2e_ups_qualifiers}:+prof:+p383b )
      setupRequired( KinKal  ${kinkal_ver} -q +${mu2e_ups_qualifiers}:+prof:+p383b )
      setupRequired( heppdt  ${heppdt_ver} -q +${mu2e_ups_qualifiers}:+prof )
      setupRequired( xerces_c  ${xerces_ver} -q +${mu2e_ups_qualifiers}:+prof )

Flavor     = ANY
Qualifiers = "${qualifiers_value}:trig:debug"

  Action = GetFQDir
      envSet( \${UPS_PROD_NAME_UC}_FS, "" )
      setupRequired( cetpkgsupport )
      execute( "get-directory-name subdir", NO_UPS_ENV, \${UPS_PROD_NAME_UC}_FS )
      envSet (OFFLINE_FQ, \${\${UPS_PROD_NAME_UC}_FS}.${offline_fq_value}.trig.debug)
      envSet (MU2E_UPS_QUALIFIERS, +${mu2e_ups_qualifiers}:+debug )
      setupRequired( art  ${art_ver} -q +${mu2e_ups_qualifiers}:+debug )
      setupRequired( BTrk  ${btrk_ver} -q +${mu2e_ups_qualifiers}:+debug:+p383b )
      setupRequired( KinKal  ${kinkal_ver} -q +${mu2e_ups_qualifiers}:+debug:+p383b )
      setupRequired( heppdt  ${heppdt_ver} -q +${mu2e_ups_qualifiers}:+debug )
      setupRequired( xerces_c  ${xerces_ver} -q +${mu2e_ups_qualifiers}:+debug )
      setupRequired( mu2e_artdaq_core ${mu2e_artdaq_core_ver} -q ${mu2e_ups_qualifiers}:+${MU2E_ART_SQUALIFIER}:+debug )

Flavor     = ANY
Qualifiers = "${qualifiers_value}:trig:prof"

  Action = GetFQDir
      envSet( \${UPS_PROD_NAME_UC}_FS, "" )
      setupRequired( cetpkgsupport )
      execute( "get-directory-name subdir", NO_UPS_ENV, \${UPS_PROD_NAME_UC}_FS )
      envSet (OFFLINE_FQ, \${\${UPS_PROD_NAME_UC}_FS}.${offline_fq_value}.trig.prof)
      envSet (MU2E_UPS_QUALIFIERS, +\${mu2e_ups_qualifiers}:+prof )
      setupRequired( art  ${art_ver} -q +${mu2e_ups_qualifiers}:+prof )
      setupRequired( BTrk  ${btrk_ver} -q +${mu2e_ups_qualifiers}:+prof:+p383b )
      setupRequired( KinKal  ${kinkal_ver} -q +${mu2e_ups_qualifiers}:+prof:+p383b )
      setupRequired( heppdt  ${heppdt_ver} -q +${mu2e_ups_qualifiers}:+prof )
      setupRequired( xerces_c  ${xerces_ver} -q +${mu2e_ups_qualifiers}:+prof )
      setupRequired( mu2e_artdaq_core ${mu2e_artdaq_core_ver} -q ${mu2e_ups_qualifiers}:+${MU2E_ART_SQUALIFIER}:+prof )

Common:
  Action = setup
    prodDir()
    setupEnv()
    envSet (OFFLINE_INC, \${OFFLINE_DIR}/include/Offline )
    envSet (MU2E_BASE_RELEASE, \${OFFLINE_DIR} )
    envSet (MU2E_SEARCH_PATH, \${OFFLINE_DIR}/config/Offline:\${MU2E_DATA_PATH} )
    envSet (FHICL_FILE_PATH, \${OFFLINE_DIR}/config/Offline:\${OFFLINE_DIR}/config/Offline/fcl )
    envSet (OFFLINE_VERSION, \${UPS_PROD_VERSION} )
    pathAppend( ROOT_INCLUDE_PATH, \${OFFLINE_INC})

    exeActionRequired(GetFQDir)
    envSet (OFFLINE_LIB, \${OFFLINE_DIR}/\${OFFLINE_FQ}/lib )
    pathAppend( CET_PLUGIN_PATH, \${OFFLINE_LIB} )

    if ( test \`uname\` = "Darwin" )
      pathPrepend(DYLD_LIBRARY_PATH, \${\${UPS_PROD_NAME_UC}_LIB})
    else()
      pathPrepend(LD_LIBRARY_PATH, \${\${UPS_PROD_NAME_UC}_LIB})
    endif ( test \`uname\` = "Darwin" )

End:
# End Group definition
#*************************************************
EOF

unset art_ver
unset heppdt_ver
unset btrk_ver
unset kinkal_ver
unset xerces_ver
unset mu2e_artdaq_core_ver
unset qualifiers_value
unset offline_fq_value
unset mu2e_ups_qualifiers
