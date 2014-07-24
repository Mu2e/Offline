#!/usr/bin/env bash
#----------------------------------------------------------------------------------
#  call : run_stnmaker.sh job_number -j job_name -n [nevents[:d]] -R caf:zee_pyt
#  example (type in CafGui window): run_stnfilter.sh $ -j etau08-strip -n 5000
#  if running on CAF WORK_DIR=CAF and no cd is required
#
#  remember that "-s" doesn't want to work
#
# call: run.sh  job_number \
#                            -v                         \
#                            -d output_dataset_name     \
#                            -D input_dataset_name      \
#                            -e exe_file                \
#                            -E full_exe_file           \
#                            -f fileset                 \
#                            -i input_tcl_file          \
#                            -J TCLFILE                 \
#                            -j job_name[:mc]           \
#                            -n nevents[:x]             \
#                            -o output_dir              \
#                            -p remote_server:password  \
#                            -w working_directory       \
#                            -x env_var=value
#
#  using jobs/tcl/stnmaker/job_name.tcl as  the main TCL file
#        input_tcl_file            to define input
#        output file name if       job_name.stntuple
#
#  using tiki: add "-x USE_TIKI=hostname"
#
#
#  parameters:
#
# -v: turn on script debug mode, should go first
#
# "x" = "t" or "d" or "p"
# t - use Totalview
# d - use GDB
#
# to run a typical stripping job:
# -------------------------------
# Cafcom stnmaker.tgz icaf:stnmaker-strip-gjet08.$.tgz murat@fnal.gov short 101 114 \
# "./bin/$BFARCH/run.sh -J $ -j stnmaker.strip-gjet08 -D gjet08 -d gjet08 \
#                         -o caf:strip-gjet08"
#
# to run a stnmaker_prod job on a dataset eexo18:
# -----------------------------------------------
# Cafcom stnmaker_prod.tgz icaf:eexo18.$.tgz murat@fnal.gov medium 1 5 \
# "./cdfopr/scripts/run.sh -J $ -j stnmaker.prod -D eexo18 -d eexo18 \
#  -e bin/$BFARCH/stnmaker_prod.exe \
#  -o caf:eexo18"
#
# to run a stnmaker_prod job on a dataset eexo08 and copy the output to fcdfdata030:
# ---------------------------------------------------------------------------------
# Cafcom stnmaker_prod.tgz icaf:eexo08.$.tgz murat@fnal.gov medium_dcache 1 5 \
# "./bin/$BFARCH/run.sh -J $ -j stnmaker.prod -i jobs/tcl/stnmaker/prereq-taumet.tcl \
#  -D eexo08 -d eexo08 -e stnmaker_prod.exe \
#  -o ewk@fcdfdata030:/cdf/scratch/ewk/eexo08"
#
# to run a stnmaker_prod job on a dataset bhel08, filter high-Pt electrons 
# and copy the output to fcdfdata030:
# ---------------------------------------------------------------------------------
# Cafcom stnmaker_prod.tgz icaf:bhel08.$.tgz murat@fnal.gov medium_dcache 1 70 \
# "./bin/$BFARCH/run.sh -J $ -j stnmaker.prod -i jobs/tcl/stnmaker/prereq-hpte.tcl \
# -D bhel08 -d bhel08 -e stnmaker_prod.exe -o ewk@fcdfdata030:/cdf/scratch/ewk/bhel08"
#
# test in interactive mode:
# --------------------------
# run.sh 75 -j stnmaker.strip-etau08 -D etau08 -d etau08 -e stnmaker.exe \
#           -o $PWD/results/stnmaker-strip-etau08 -n 100
#------------------------------------------------------------------------------
. Stntuple/scripts/common_procedures

echo [run.sh:$LINENO]: Stntuple/scripts/parse_parameters $*
. Stntuple/scripts/parse_parameters $*
export          OFFVER=`cat $WORK_DIR/.base_release`
export         OFFLINE=`cat $WORK_DIR/.base_release | sed s'/\.//g'` 
#------------------------------------------------------------------------------
#  so far - a kludge: if SOURCE_ME is defined, source $WORK_DIR/source_me
#  need to load proper shared libraries for stnfit.exe
#  if CDFSOFT env var can be defined, use it, otherwise assume that CDF 
#  software is not available and source local 'source_me' file
#------------------------------------------------------------------------------
echo [run.sh:$LINENO] \$SOURCE_ME=$SOURCE_ME

if [ -f $WORK_DIR/source_me ] ; then 
  source $WORK_DIR/source_me 
else 
  if [ ".`echo $CDFSOFT`" != '.' ] ; then
    echo MU2E skip     source $CDFSOFT/cdf2.shrc
    echo [run.sh:$LINENO] setup cdfsoft2 $OFFVER
    echo MU2E skip setup cdfsoft2 $OFFVER
  else
    source cdfopr/scripts/source_me.sh
  fi
fi

if [ .$USE_TIKI != "." ] ; then
  echo [run.sh]: ... using TIKI  JOB_NAME=$JOB_NAME
#-----------------------------------------------------------------------
# read parameters from the tiki database, remember that in this case
# input file is also stored in tiki
#-----------------------------------------------------------------------
#  tiki_host=$USE_TIKI
#  parameter_file=$DATASET.$FULL_JOB_NAME.parameters
#  delimitor=concat.parameters

#  echo [run.sh.tiki]: parameter_file=$parameter_file

#  echo [run.sh.tiki]: cdfopr/scripts/tiki_get_data $tiki_host wiki_page . $BOOK.$DATASET \
#                      $delimitor parameter_file=$parameter_file
#  . cdfopr/scripts/tiki_get_data $tiki_host wiki_page . $BOOK.$DATASET $delimitor > $parameter_file
#  cat $parameter_file

#  script=cdfopr/scripts/tiki_get_parameter
#  input_tcl_file=`$script $parameter_file input_tcl_file $JOB_NUMBER`
#     output_node=`$script $parameter_file output_node    $JOB_NUMBER`
#      output_dir=`$script $parameter_file output_dir     $JOB_NUMBER`

#  echo [run.sh.tiki]: input_tcl_file=$input_tcl_file
#  echo [run.sh.tiki]:    output_node=$output_node
#  echo [run.sh.tiki]:     output_dir=$output_dir

#  cmd=". cdfopr/scripts/get_file tiki:$tiki_host $WORK_DIR $BOOK $DATASET $input_tcl_file ${WORK_DIR}/$input_tcl_file"
#  echo [run.sh.tiki]: $cmd ;

#  $cmd

#  export INPUT_TCL_FILE=${WORK_DIR}/$input_tcl_file

#  echo [run.sh.tiki]: done, INPUT_TCL_FILE=$INPUT_TCL_FILE 
#  . cdfopr/scripts/parse_output_dir ${output_node}:$output_dir  $DEBUG_SCRIPT
#  export RUSER=$OP1
#  export RHOST=$OP2
#  export  RDIR=$OP3
#  echo [run.sh.tiki]: RUSER.RHOST.RDIR.tiki. = $RUSER.$RHOST.$RDIR
fi


echo [run.sh:$LINENO]: JOB_OUTPUT_DIR=$JOB_OUTPUT_DIR
echo [run.sh:$LINENO]: PATH=$PATH
echo [run.sh:$LINENO]: which rcp = `which rcp`
#------------------------------------------------------------------------------
export JOB_TCL_DIR=$WORK_DIR/$JOB_TCL_DIR

index=`printf "%04i" $JOB_NUMBER`

exe_stub=`echo $EXE | awk -F . '{print $1}'`

if [ .$JOB_NAME != "." ] ; then 
  if [ ".$SCRIPT" == "." ] ; then
    name_stub=${exe_stub}-${JOB_NAME}
  else
    name_stub=${exe_stub}
  fi
fi

export DSID=`echo $DATASET | sed 's#/#-#g'`
if [ .$TIKI_DSID == "." ] ; then
  export TIKI_DSID=$DSID
fi

name_stub=${name_stub}-$TIKI_DSID

if [ .$JOB_OUTPUT_DIR == "." ] ; then
  export JOB_OUTPUT_DIR=${WORK_DIR}/results/${name_stub}.${index}
fi

if [ $DEBUG_SCRIPT != 0 ] ; then 
  echo [run.sh:$LINENO]: exe_stub=$exe_stub
  echo [run.sh:$LINENO]: JOB_NAME=$JOB_NAME
  echo [run.sh:$LINENO]: name_stub=$exe_stub
  echo [run.sh:$LINENO]: JOB_OUTPUT_DIR=$JOB_OUTPUT_DIR
fi


if [ .$TCLFILE == "." ] ; then
  tcl_stub=$exe_stub
  if [ .$JOB_TCL != "." ] ; then 
    tcl_stub=${JOB_TCL}
  fi
 
  if [ .$JOB_TCL == "." ] ; then
    export TCLFILE=${JOB_TCL_DIR}/${exe_stub}/${tcl_stub}.tcl
  else
    export TCLFILE=${JOB_TCL_DIR}/${TCL_SUBDIR}/${tcl_stub}.tcl
  fi
fi

cd  $WORK_DIR

if [ ! -f $WORK_DIR/source_me ] ; then
  if [ ".$SRT_QUAL" == "." ] ; then SRT_QUAL=default ; fi
  echo [run.sh:$LINENO]: MU2E skip srt_setup -a SRT_QUAL=$SRT_QUAL
  echo [run.sh:$LINENO]: MU2E skip setup fcp
fi

export LD_LIBRARY_PATH=/usr/krb5/lib:$LD_LIBRARY_PATH
# echo LD_LIBRARY_PATH=$LD_LIBRARY_PATH

echo [run.sh:$LINENO]: SKIP export PRODUCTS=/cdf/home/cdfopr/upsdb:$PRODUCTS

if [ .$OUTPUT_FILE != "." ] ; then nm=${OUTPUT_FILE}.${index}
                              else nm=${name_stub}.${index}
fi
#-----------------------------------------------------------------------
# avoid name clashes in case both - output file and stntuple are defined
# (normally this should not be happening, just in case)
#-----------------------------------------------------------------------
export     OUTPUT_FILE=${JOB_OUTPUT_DIR}/${nm}.s
export OUTPUT_STNTUPLE=${JOB_OUTPUT_DIR}/${nm}.s
export         LOGFILE=${JOB_OUTPUT_DIR}/${nm}.log

export OUTPUT_NAME=$nm

if [ .$EXEFILE == "." ] ; then
  export       EXEFILE=${WORK_DIR}/bin/$BFARCH/${EXE}
fi

if [ ! -f $JOB_OUTPUT_DIR ] ; then mkdir -p $JOB_OUTPUT_DIR ; fi
cd $JOB_OUTPUT_DIR

if [ $DEBUG_SCRIPT != 0 ] ; then 
  echo [run.sh:$LINENO]: OUTPUT_STNTUPLE = $OUTPUT_STNTUPLE
  echo [run.sh:$LINENO]: OUTPUT_FILE     = $OUTPUT_FILE
  echo [run.sh:$LINENO]: LOGFILE         = $LOGFILE
  echo [run.sh:$LINENO]: EXEFILE         = $EXEFILE
  echo [run.sh:$LINENO]: TCLFILE         = $TCLFILE
  echo [run.sh:$LINENO]: SCRIPT          = $SCRIPT
  echo [run.sh:$LINENO]: PARAMETERS      = $PARAMETERS
fi
#-----------------------------------------------------------------------
#  initialize the logfile
#-----------------------------------------------------------------------
. $WORK_DIR/Stntuple/scripts/init_logfile


# if [ ! -f libreadline.so.4.1 ] ; then 
#   ln -s /usr/lib/libreadline.so $JOB_OUTPUT_DIR/libreadline.so.4.1 
# fi
# 
# if [ ! -f libhistory.so.4.1 ] ; then 
#   ln -s /usr/lib/libhistory.so $JOB_OUTPUT_DIR/libhistory.so.4.1 
# fi

export LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH

.  $WORK_DIR/Stntuple/scripts/submit_job
echo [run.sh]: "RETURN_CODE=$RETURN_CODE"                                    >> $LOGFILE
#------------------------------------------------------------------------------
#  after job finished print some statistics
#------------------------------------------------------------------------------
echo [run.sh:$LINENO]: "###########################################################" >> $LOGFILE
echo [run.sh:$LINENO]: "# HOST: `hostname -i`" " `hostname -f`"                      >> $LOGFILE
echo [run.sh:$LINENO]: "# " `cat /proc/cpuinfo | \
                     egrep processor\|model\|cpu\|cache | uniq`                      >> $LOGFILE
echo [run.sh:$LINENO]: "###########################################################" >> $LOGFILE
echo [run.sh:$LINENO]: "# my process ID: " $$                                        >> $LOGFILE
echo [run.sh:$LINENO]: "# " `cat /proc/$$/stat`                                      >> $LOGFILE
echo [run.sh:$LINENO]: "# CPU time (usr/sys, sec): " \
      `cat /proc/$$/stat | awk '{print $16/100 "\/" $17/100}'`                       >> $LOGFILE
echo [run.sh:$LINENO]: "# starting cleanup actions on $HOSTNAME                    " >> $LOGFILE
echo [run.sh:$LINENO]: "###########################################################" >> $LOGFILE
ls  -alF $JOB_OUTPUT_DIR                                                             >> $LOGFILE
echo [run.sh:$LINENO]: "###########################################################" >> $LOGFILE
#------------------------------------------------------------------------------
#  final actions - clean up the WORK_DIR if running on CAF
#------------------------------------------------------------------------------
echo [run.sh:$LINENO]: MU2E skip . $WORK_DIR/cdfopr/scripts/do_cleanup

