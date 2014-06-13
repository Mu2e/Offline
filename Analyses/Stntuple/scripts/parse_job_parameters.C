//-----------------------------------------------------------------------------
// initialize default values of the job parameters
// parse command line (Parameters) to redefine them 
// recognized keys: 
//  /job_name= 
//  /job_number= 
//  /l3_trig_path=
//  /mc=
//  /calib_pass=
//  /debug=
//  /little=
//  /grl=                        (good run list)
//  /output=                     (name of the output file)
//-----------------------------------------------------------------------------
#include "global_vars.h"

int parse_job_parameters(TString& Parameters, StnAnaGlobals_t& Glob) {
  TString  fileset, job_number;
  int      loc, loc1, loc2, len;

  printf("job = %s\n",Parameters.Data());
//-----------------------------------------------------------------------------
// environment variables which may be defined by default
// tau_ana keeps list of the known jobs
//-----------------------------------------------------------------------------
  fileset      = gSystem->Getenv("FILESET");
  g.JobName    = gSystem->Getenv("JOB_NAME");
  g.JobNumber  = gSystem->Getenv("JOB_NUMBER");
  g.CalibPass  = gSystem->Getenv("CALIB_PASS");
  g.L3TrigPath = gSystem->Getenv("L3_TRIG_PATH");
//-----------------------------------------------------------------------------
// job_name: name of the job (may be redefined)
//-----------------------------------------------------------------------------
  loc     = Parameters.Index("/job_name=");
  len     = Parameters.Length();
  if (loc > 0)  {
    loc1 = loc+strlen   ("/job_name=");
    loc2 = Parameters.Index("/",loc+1);
    if (loc2 < 0) loc2 = len;
    g.JobName = Parameters(loc1,loc2-loc1);

    printf("JOB name    : %s\n",g.JobName.Data());
  }

  loc = Parameters.Index("/");
  if (loc > 0)  g.Task = Parameters(0,loc);
  else          g.Task = Parameters;

  if (g.JobName == "") {
    loc = Parameters.Index("(");
    g.JobName = Parameters(0,loc);
  }
//-----------------------------------------------------------------------------
// job_number:
//-----------------------------------------------------------------------------
  loc     = Parameters.Index("/job_number=");
  if (loc > 0)  {
    loc1 = loc+strlen   ("/job_number=");
    loc2 = Parameters.Index("/",loc+1);
    if (loc2 < 0) loc2 = len;
    g.JobNumber = Parameters(loc1,loc2-loc1);
  }
  gSystem->Setenv("JOB_NUMBER",g.JobNumber.Data());
  printf("%-15s : %s\n","JOB_NUMBER",g.JobNumber.Data());
//-----------------------------------------------------------------------------
// l3_trig_path:
// g.L3TrigPath is processed by murat/ana/scripts/setup_trigger_path.C called
//              from stnana.C
// I didn't find which part of the code (if any) deals with 
// the environment variable L3_TRIG_PATH
//-----------------------------------------------------------------------------
  loc     = Parameters.Index("/l3_trig_path=");
  if (loc > 0)  {
    loc1 = loc+strlen   ("/l3_trig_path=");
    loc2 = Parameters.Index("/",loc+1);
    if (loc2 < 0) loc2 = len;
    g.L3TrigPath = Parameters(loc1,loc2-loc1);
  }
  printf("%-15s : %s\n","L3_TRIG_PATH",g.L3TrigPath.Data());
//-----------------------------------------------------------------------------
// MC flag (default: -1) this flag set to 1 indicates that for a MC dataset 
// some MC-specific actions need to be taken
//-----------------------------------------------------------------------------
  loc = Parameters.Index("/mc");
  if (loc > 0)  {
    loc = Parameters.Index("/mc=");
    if (loc <= 0) g.DoMc = 1;
    else {
      loc1 = loc+strlen   ("/mc=");
      loc2 = Parameters.Index("/",loc+1);
      if (loc2 < 0) loc2 = len;
      g.DoMc = atoi(Parameters(loc1,loc2-loc1).Data());
    }
  }
//-----------------------------------------------------------------------------
// calibration pass - set reasonable defaults, which can be explicitly redefined
// if calib_pass has not been specified via environment, set defaults.
// default pass is '05' for the data and '05.mc' for MC
//-----------------------------------------------------------------------------
  if (g.CalibPass == "") {
    if (g.DoMc > 0) g.CalibPass = "05.mc";
    else            g.CalibPass = gEnv->GetValue("Stnana.CalibPass","05");
  }
				        // allow calib pass to be redefined on 
					// the command line - check for it
					// preserve backward compatibility:
					// '/pass'
  loc = Parameters.Index("/calib_pass=");
  if (loc > 0) {
    loc1 = loc+strlen   ("/calib_pass=");
    loc2 = Parameters.Index("/",loc+1);
    if (loc2 < 0) loc2 = len;
    g.CalibPass = Parameters(loc1,loc2-loc1);
  }
  else {
    loc = Parameters.Index("/pass=");
    if (loc > 0)  {
      loc1 = loc+strlen   ("/pass=");
      loc2 = Parameters.Index("/",loc+1);
      if (loc2 < 0) loc2 = len;
      g.CalibPass = Parameters(loc1,loc2-loc1);
    }
  }
//-----------------------------------------------------------------------------
// debug
//-----------------------------------------------------------------------------
  loc         = Parameters.Index("/debug");
  if (loc > 0) g.Debug = 1;
  else         g.Debug = 0;
//-----------------------------------------------------------------------------
// flag controlling writing out "little" ntuple
//-----------------------------------------------------------------------------
  loc     = Parameters.Index("/little=");
  if (loc > 0)  {
    g.DoLittle = 1;
    loc1 = loc+strlen   ("/little=");
    loc2 = Parameters.Index("/",loc+1);
    if (loc2 < 0) loc2 = len;
    g.LittleFileName = Parameters(loc1,loc2-loc1);

    printf("%-15s : %s\n","LITTLE_FILE",g.LittleFileName.Data());
  }
//-----------------------------------------------------------------------------
// flag controlling writing output of stripping
//-----------------------------------------------------------------------------
  g.OutputFileName = "";
  loc = Parameters.Index("/output");
  if (loc > 0)  {
    loc1 = loc+strlen("/output=");
    loc2 = Parameters.Index("/",loc+1);
    if (loc2 < 0) loc2 = len;
    if (loc2 > loc1) {
      g.OutputFileName = Parameters(loc1,loc2-loc1);
    }
    else {
      g.OutputFileName = Form("%s",gSystem->Getenv("DSID"));

      if (fileset  != "") g.OutputFileName += Form(".%s",fileset.Data());
      g.OutputFileName += Form(".%s",g.JobName.Data());
      g.OutputFileName += ".stn";
    }
    printf("%-15s : %s\n","OUTPUT_FILE",g.OutputFileName.Data());
  }
//-----------------------------------------------------------------------------
// saving the histograms: /save
//-----------------------------------------------------------------------------
  g.HistFileName = "";
  loc = Parameters.Index("/save");
  if (loc > 0)  {
    loc1 = loc+strlen   ("/save=");
    loc2 = Parameters.Index("/",loc+1);
    if (loc2 <= 0) loc2 = len;
    if (loc2 > loc1) {
      g.HistFileName = Parameters(loc1,loc2-loc1);
    }
    else {
      g.HistFileName   = Form("%s",gSystem->Getenv("DSID"));

      if (fileset  != "") g.HistFileName += Form(".%s",fileset.Data());
      g.HistFileName += Form(".%s",g.JobName.Data());
      g.HistFileName += ".hist";
    }

    printf("%-15s : %s\n","HISTOGRAM_FILE",g.HistFileName.Data());
  }
//-----------------------------------------------------------------------------
// handle good run list  ( use -x STNTUPLE_GRL=ETF,141544,156487 )
// use GRL only for non-MC jobs
// finally allow to redefine GRL from the command line
//-----------------------------------------------------------------------------
  const char    *env, good_run_list[200];
  TString       s("none,1,999999999");
  env = gSystem->Getenv("STNTUPLE_GRL");
  if (env != 0) s = env;

  loc = Parameters.Index("/grl=");
  //  printf("loc= %i\n",loc);
  if (loc > 0) {
    loc1 = loc+strlen   ("/grl=");
    loc2 = Parameters.Index("/",loc+1);
    if (loc2 <= 0) loc2 = len;
    //    printf("loc1,loc2,len= %i %i %i\n",loc1,loc2,len);
    s = Parameters(loc1,loc2-loc1);
    //    printf("s= %s\n",s.Data());
  }

  s.ReplaceAll(","," ");
  sscanf(s.Data(),"%s %i %i",good_run_list,&g.MinRun,&g.MaxRun);
  
  g.GoodRunList = good_run_list;
  g.GoodRunList.ToUpper();
  printf("GOOD RUN LIST : \"%s\" with run range %i-%i\n"
	 ,g.GoodRunList.Data(),g.MinRun,g.MaxRun);
//-----------------------------------------------------------------------------
// Steve's "newcuts" flag
//-----------------------------------------------------------------------------
  loc = Parameters.Index("/newcuts");
  if (loc > 0) g.IDMode = 2;   // use the new cuts
  else         g.IDMode = 1;   // use the old cuts

  return 0;
}
