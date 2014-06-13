///////////////////////////////////////////////////////////////////////////////
//  STNTUPLE implementeations of various good run lists
//-----------------------------------------------------------------------------
//  calculate lumi for the ETF definition of a good run list
//  integrated luminosity:     all runs          good runs
//  ---------------------------------------------------------------------------
//  runs 141544-154799:       85630.539062     50814.324219
//  runs 141544-156487:      106742.226562     70732.757812
//
// 2009-05-23: for partially good runs consider only their sections when the 
//             whole detector has been operational
//-----------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////
#include "obj/TStnRunSummary.hh"
#include "obj/TStnGoodRunList.hh"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TEnv.h"
#include "TUrl.h"
//-----------------------------------------------------------------------------
//  variables local for the file
//-----------------------------------------------------------------------------
namespace {
  TStnGoodRunList* fgGoodRunList = NULL;
}

ClassImp(TStnGoodRunList)
//-----------------------------------------------------------------------------
TStnGoodRunList::TStnGoodRunList(const char* Name, const char* Filename) : 
  fMinGoodRunLum(1.)
  , fNRanges(0)
  , fUseUncheckedRuns(1)
{
  TString name(Name);

  name.ToUpper();

  SetName(name.Data());
  SetTitle(name.Data());

  fRunSummary   = new TStnRunSummary();
  fFile         = NULL;
  fTree         = NULL;
  fMinRunNumber = -1;
  fMaxRunNumber = -1;
  fNEntries     = -1;
  fCurrentEntry = -1;

  Init(Filename);

  fgGoodRunList = this;

  if      (fName == "ETF"   ) {
    SetGoodRunRoutine(GoodRunListEtf    );
    printf(" **** ETF good run list is used\n");
  }
  else if (fName == "WZ_PRD") {
    SetGoodRunRoutine(GoodRunList_WZ_PRD);
    printf(" **** WZ_PRD good run list is used\n");
  }
  else if (fName == "DQM_V6") {
    SetGoodRunRoutine(GoodRunList_DQM_V6);
    printf(" **** DQM V6 good run list is used\n");
  }
  else if (fName.Index("DQM_V7") == 0) {
//-----------------------------------------------------------------------------
//  DQM V7 lists: DQM_V7:XYZ
//  X: electron bit:  0: no requirements, 1: require CES/PES
//  Y: muon     bit:  0: no requirements
//                    1: require CMUP+CMX (default)
//                    2: require CMUP, ignore CMX
//                    3: require CMX , ignore CMUP
//  Z: silicon  bit:  0: ignore silicon
//                    1: require silicon (SVX+ISL+L00)
//-----------------------------------------------------------------------------
    SetGoodRunRoutine(GoodRunList_DQM_V7);
    printf(" **** DQM V7 good run list is used\n");
  }
  else if (fName.Index("DQM_V13") == 0) {
//-----------------------------------------------------------------------------
//  DQM V13 lists: DQM_V13:XYZ , same as DQM_V7
//-----------------------------------------------------------------------------
    SetGoodRunRoutine(GoodRunList_DQM_V13);
    printf(" **** DQM V13 good run list is used\n");
  }
  else if (fName.Index("DQM_V27") == 0) {
    SetGoodRunRoutine(GoodRunList_DQM_V27);
    printf(" **** DQM V27 good run list from Apr 14 2009 is used\n");
  }
  else if (fName.Index("DQM_V32") == 0) {
    SetGoodRunRoutine(GoodRunList_DQM_V32);
    printf(" **** DQM V32 good run list from Mar 12 2010 is used\n");
  }
  else if (fName.Index("DQM_V34") == 0) {
    SetGoodRunRoutine(GoodRunList_DQM_V34);
    printf(" **** DQM V34 good run list from May 14 2010 is used\n");
  }
  else if (fName.Index("MC_1001") == 0) {
    SetGoodRunRoutine(GoodRunList_MC_1001);
    printf(" **** MC_1001 good run list is used\n");
  }
  else {
    SetGoodRunRoutine(NULL);
    printf(" **** no good run list is used\n");
  }

  TString s1;
  int loc = fName.Index(':');

  if (loc >= 0) s1 = fName(loc+1,3);
  else          s1 = "111";

  fElectronFlag = int(s1(0)) - int('0');
  fMuonFlag     = int(s1(1)) - int('0');
  fSiliconFlag  = int(s1(2)) - int('0');
}

//-----------------------------------------------------------------------------
TStnGoodRunList::~TStnGoodRunList() {
  delete fRunSummary;
  if (fFile) delete fFile;
  fgGoodRunList = 0;
}

//-----------------------------------------------------------------------------
void TStnGoodRunList::SetListOfRuns(Int_t* List) {
//-----------------------------------------------------------------------------
// count number of runs, n is supposed to be an even number,
// List[2*n] and List[2*n+1] give lover and upper boundaries of the i-th
// range of runs
//-----------------------------------------------------------------------------
  int n=0;
  while (List[n] != -1) n++;
  fListOfRuns.Set(n,List);
  fNRanges = n/2;
}

//-----------------------------------------------------------------------------
TStnRunSummary* TStnGoodRunList::GetEntry(Int_t IEntry) {
  fTree->GetEntry(IEntry);
  return fRunSummary;
}

//-----------------------------------------------------------------------------
TStnRunSummary* TStnGoodRunList::GetRunSummary(Int_t RunNumber) {

  int             rn, rn_min, rn_max;
  double          i, isave, imin, imax; 
  TStnRunSummary* rs = 0;

  imin       = 0; 
  imax       = fNEntries-1;
  rn_min     = fMinRunNumber;
  rn_max     = fMaxRunNumber;
  i          = fCurrentEntry;
  rn         = fRunSummary->RunNumber();


//    printf(" RunNumber         = %8i\n",RunNumber);
//    printf(" gMinRunNumber     = %8i\n",gMinRunNumber);
//    printf(" gMaxRunNumber     = %8i\n",gMaxRunNumber);
//    printf(" gMaxCurrentNumber = %8i\n",gCurrentRunNumber);

  if ((RunNumber < rn_min) || (RunNumber > rn_max)) {
    Error("GetRunSummary",Form("Run %8i outside the range",RunNumber));
    return NULL;
  }
  else if (RunNumber == rn_min) {
    fTree->GetEntry((int) imin);
    fCurrentEntry = imin;
    rs            = fRunSummary;
    rn            = RunNumber;
    goto END;
  }
  else if (RunNumber == rn_max) {
    fTree->GetEntry((int) imax);
    fCurrentEntry = imax;
    rs            = fRunSummary;
    rn            = RunNumber;
    goto END;
  }
  else if (RunNumber == rn) {
    rs = fRunSummary;
    goto END;
  }
//-----------------------------------------------------------------------------
// requested run number is somewhere in between, not the first and not the
// last records
//-----------------------------------------------------------------------------
  if (RunNumber < rn) {
    imax   = i;
    rn_max = rn;
  }
  else {
    imin   = i;
    rn_min = rn;
  }

  isave = imin;
  i     = (imin+imax+1)/2;

  while (i != isave) {
    isave = i;
    
    fTree->GetEntry((int) i);
    rn            = fRunSummary->RunNumber();
    fCurrentEntry = i;

//      printf(" --- i, imin, imax, rn_min, rn_max, run = %6i %6i %6i %8i %8i %8i\n",
//  	   i,imin, imax,rn_min, rn_max,rn);

    if (rn == RunNumber) {
      rs = fRunSummary;
      break;
    }
    else if (imax == imin+1) {
      imax   = imin;
      rn_max = rn_min;
    }
    else if (rn > RunNumber) {
      imax   = i;
      rn_max = rn;
    }
    else {
      imin   = i;
      rn_min = rn;
    }
    i     = (imin+imax+1)/2;
  }

 END:;

  //  printf(" run , found = %8i %3i\n",rn, (rs != 0));
  return rs;
}


//-----------------------------------------------------------------------------
Int_t TStnGoodRunList::Init(const char* Filename) {
  int rc;

  TString     cmd;
  char        filename[200];
  const char* fn;


  if (Filename != 0) {
//-----------------------------------------------------------------------------
// filename has been specified explicitly
//-----------------------------------------------------------------------------
    strcpy(filename,Filename);
  }
  else {
//-----------------------------------------------------------------------------
//  defaults ... start from looking up in .rootrc files, use the last one
//-----------------------------------------------------------------------------
    fn = gEnv->GetValue("Stntuple.RunSummary","undefined");

    if (strcmp(fn,"undefined") != 0) {
//-----------------------------------------------------------------------------
// assume web-based access, do WGET 
//-----------------------------------------------------------------------------
// ##
// #

//       cmd = Form("/usr/bin/wget -O `echo %s | awq '{%s %s", 
// 		 fn.Data(),book,dset,fset->GetName());

//       TString grep_filesets;

//       if (fPrintLevel > 0)
// 	printf(" DEBUG: Retrieving data with: %s\n",cmd.Data());
//       pipe = gSystem->OpenPipe(cmd,"r");
//       int first = 1;
//       while (fgets(buf,10000,pipe)) { 
	

// ###
      strcpy(filename,fn);
    }
    else {
//-----------------------------------------------------------------------------
//  nothing in .rootrc files, check environment
//-----------------------------------------------------------------------------
      fn = gSystem->Getenv("STNTUPLE_RUN_SUMMARY");
//-----------------------------------------------------------------------------
//  finally fall back to the default one (stored on FCDFLNX3)
//-----------------------------------------------------------------------------
      if (fn != 0) {
	strcpy(filename,fn);
      }
      else {
	sprintf(filename,"root://fcdfdata122.fnal.gov/%s",
		"/export/data2/ewk/run_summary/run_summary.root");
      }
    }
  }
//-----------------------------------------------------------------------------
//  at this point filename is guaranteed to be defined
//-----------------------------------------------------------------------------
  TUrl    url(filename);
  TString dcache_io_tunnel = gSystem->Getenv("DCACHE_IO_TUNNEL");

  if (strcmp(url.GetProtocol(),"dcap") == 0) {
//-----------------------------------------------------------------------------
// assume this is diskpool, force unsetenv 
//-----------------------------------------------------------------------------
    gSystem->Unsetenv("DCACHE_IO_TUNNEL");  
  }

  fFile = TFile::Open(filename);

  if (! fFile->IsOpen()) {
    Error("Init",Form("can\'t open file %s",filename));
    delete fFile;
    rc = -1;
  }
  else {
    fTree = (TTree*) fFile->Get("run_summary");
    fTree->SetBranchAddress("RunSummary",&fRunSummary);

    fTree->GetEntry(0);
    fMinRunNumber = fRunSummary->RunNumber();
    fNEntries     = fTree->GetEntries();

    fCurrentEntry = fNEntries-1;
    fTree->GetEntry((int) fCurrentEntry);
    fMaxRunNumber = fRunSummary->RunNumber();

    printf(" TStnGoodRunList::Init: use %s to initialize the good run list\n",
	   filename);
    rc = 0;
  }

  gSystem->Setenv("DCACHE_IO_TUNNEL",dcache_io_tunnel.Data());

  return rc;
}

//_____________________________________________________________________________
Int_t TStnGoodRunList::GoodRunListEtf(int RunNumber, int RunSection, int Mask) {
  // today, Jan 25 2003 ETF good run list is defined only for runs <= 154799
  // call with RunNumber < 0 sets up print level

  TStnGoodRunList* grl = NULL;

  int electron_flag, muon_flag, good_run(0); // , silicon_flag

  TString*        tt;
  TStnRunSummary* rs;

//-----------------------------------------------------------------------------
// initialization
//-----------------------------------------------------------------------------
  grl = fgGoodRunList;
  
  if (RunNumber < 141544) return good_run ;

  rs = grl->GetRunSummary(RunNumber);
  if (rs != 0) {

    good_run = 1;

    good_run *= rs->ClcStatusBit();
    good_run *= rs->L1tStatusBit();
    good_run *= rs->L2tStatusBit();
    good_run *= rs->L3tStatusBit();
    good_run *= rs->CalStatusBit();

    electron_flag = grl->ElectronFlag();
    muon_flag     = grl->MuonFlag    ();
    //    silicon_flag  = grl->SiliconFlag ();

    if (electron_flag == 1) good_run *= rs->SmxStatusBit();

    good_run *= rs->CalOfflineBit();
    good_run *= rs->CotOfflineBit();
//-----------------------------------------------------------------------------
//  CMUP: online and offline bits
//-----------------------------------------------------------------------------
    if ((muon_flag == 1) || (muon_flag == 2)) {
      good_run *= rs->CmuStatusBit();
      good_run *= rs->CmpStatusBit();
      good_run *= rs->CmuOfflineBit();
      good_run *= rs->CmpOfflineBit();
    }
//-----------------------------------------------------------------------------
//  CMX : ON for runs after 150145
//-----------------------------------------------------------------------------
    if ((muon_flag == 1) || (muon_flag == 3)) {
      if (RunNumber < 150145) {
	good_run = 0;
      }
      else {
	good_run *= rs->CmxStatusBit();
	good_run *= rs->CmxOfflineBit();
      }
    }
//-----------------------------------------------------------------------------
// this is a list of already fixed exceptions
//-----------------------------------------------------------------------------
    if (RunNumber == 141989) good_run = 0;
    if (RunNumber == 142130) good_run = 1;
    if (RunNumber == 142131) good_run = 1;
    if (RunNumber == 144327) good_run = 1;
    if (RunNumber == 144680) good_run = 1;
    if (RunNumber == 144691) good_run = 1;
    if (RunNumber == 144880) good_run = 0; // COT web page: BAD, COT offl bit=1, remove
    if (RunNumber == 145613) good_run = 1; // all onl  bits null, otherwise OK, keep
    if (RunNumber == 145669) good_run = 1; // all onl  bits null, otherwise OK, keep
    if (RunNumber == 150820) good_run = 1; // CMX onl  bit 0: miniskirt, keep
    if (RunNumber == 150821) good_run = 1; // CMX onl  bit 0: miniskirt, keep
    if (RunNumber == 150856) good_run = 0; // CMX offl bit set to 1 by mistake, remove
    if (RunNumber == 151556) good_run = 1; // CMU online bit set to 0 by mistake, keep
    if (RunNumber == 151845) good_run = 1; // CAL and CMX online bits NULL, no reason, keep
    if (RunNumber == 151872) good_run = 1; // CAL onl left 0 - ptntl isolist pb, keep
    if (RunNumber == 151873) good_run = 1; // same
    if (RunNumber == 151906) good_run = 1; // same
    if (RunNumber == 151907) good_run = 1; // same
    if (RunNumber == 152967) good_run = 0; // CMX offl set to 1 by mistake, remove
    if (RunNumber == 153946) good_run = 0; // SMX HV prblms, although SMX online=1, remove
    if (RunNumber == 153985) good_run = 0; // same
//-----------------------------------------------------------------------------
//  3 exceptions from the 2nd part - ask Eric!
//-----------------------------------------------------------------------------
    if (RunNumber == 155324) good_run = 1; // Michael has not looked at it yet, keep
    if (RunNumber == 156370) good_run = 0; // special readout list trig tests, remove
    if (RunNumber == 156371) good_run = 0; // same as above
//-----------------------------------------------------------------------------
//  pre Jan'2003 shutdown runs taken with the PHYSICS_1_04 trigger table
//  (processed around Jan 20 2003, so we can start use them again)
//-----------------------------------------------------------------------------
//      if (RunNumber == 156372) good_run = 0;
//      if (RunNumber == 156401) good_run = 0;
//      if (RunNumber == 156452) good_run = 0;
//      if (RunNumber == 156457) good_run = 0;
//      if (RunNumber == 156460) good_run = 0;
//      if (RunNumber == 156464) good_run = 0;
//-----------------------------------------------------------------------------
//  Jan 13 2003: here the winter'2003 shutdown starts (last physics run 156487)
//-----------------------------------------------------------------------------
    if (RunNumber == 164844) good_run = 0;             // CSL problems
    if (RunNumber == 164870) good_run = 0;             // CSL problems
    if (RunNumber == 164871) good_run = 0;             // CSL problems
    if (RunNumber == 164872) good_run = 0;             // CSL problems
//-----------------------------------------------------------------------------
//  Sep 07 2003: fall'2003 shutdown starts (last physics run: 168892)
//  finally check the trigger table name - make sure it is PHYSICS but not TEST
//-----------------------------------------------------------------------------
    //  CHECK_TRIGGER_TABLE:
    tt = &rs->TriggerTableName();
    if ((tt->Index("PHYSICS") < 0) || (tt->Index("TEST") >= 0)) {
      good_run = 0;
    }
  }
  return good_run;
}


//_____________________________________________________________________________
Int_t TStnGoodRunList::GoodRunList_WZ_PRD(int RunNumber, int RunSection, int Mask) {
  // frozen, runs 141544-156487, do not accept "promotions"
  // total integrated luminosity = 72 pb^-1

  TStnGoodRunList* grl = NULL;

  TString*        tt;
  TStnRunSummary* rs;

  int good_run = 0;
  				// initialization
  grl = fgGoodRunList;
  
  if ((RunNumber < 141544) || (RunNumber > 156487)) return good_run ;

  rs = grl->GetRunSummary(RunNumber);
  if (rs != 0) {
//-----------------------------------------------------------------------------
// proceed with selecting the good runs: Eric's logic
//-----------------------------------------------------------------------------
    good_run = 1;

    good_run *= rs->ClcStatusBit();
    good_run *= rs->L1tStatusBit();
    good_run *= rs->L2tStatusBit();
    good_run *= rs->L3tStatusBit();
    good_run *= rs->CalStatusBit();
    good_run *= rs->CmuStatusBit();
    good_run *= rs->CmpStatusBit();
    good_run *= rs->SmxStatusBit();
//-----------------------------------------------------------------------------
//  Jun 22 2003: for post-winter'2003 shutdown runs check only online bits
//  plus COT and CAL, other bits to follow
//-----------------------------------------------------------------------------
    if (RunNumber > 165523) goto CHECK_TRIGGER_TABLE;

    good_run *= rs->CalOfflineBit();
    good_run *= rs->CotOfflineBit();

    //    if (RunNumber > 156487) goto CHECK_TRIGGER_TABLE;
//-----------------------------------------------------------------------------
//  for pre-winter'2003 shutdown runs: Eric's logic
//-----------------------------------------------------------------------------
    good_run *= rs->CmuOfflineBit();
    good_run *= rs->CmpOfflineBit();

    if (RunNumber >=150145) {
      good_run *= rs->CmxStatusBit();
      good_run *= rs->CmxOfflineBit();
    }
//-----------------------------------------------------------------------------
// this is a list of already fixed exceptions
//-----------------------------------------------------------------------------
    if (RunNumber == 141989) good_run = 0;
    if (RunNumber == 142130) good_run = 1;
    if (RunNumber == 142131) good_run = 1;
    if (RunNumber == 144327) good_run = 1;
    if (RunNumber == 144680) good_run = 1;
    if (RunNumber == 144691) good_run = 1;
    if (RunNumber == 144880) good_run = 0; // COT web page: BAD, COT offline bit=1, remove
    if (RunNumber == 145613) good_run = 1; // all online bits null, otherwise OK, keep
    if (RunNumber == 145669) good_run = 1; // all online bits null, otherwise OK, keep
    if (RunNumber == 150820) good_run = 1; // CMX online bit 0: miniskirt, keep
    if (RunNumber == 150821) good_run = 1; // CMX online bit 0: miniskirt, keep
    if (RunNumber == 150856) good_run = 0; // CMX offline bit set to 1 by mistake, remove
    if (RunNumber == 151556) good_run = 1; // CMU online bit set to 0 by mistake, keep
    if (RunNumber == 151845) good_run = 1; // CAL and CMX online bits NULL, no reason, keep
    if (RunNumber == 151872) good_run = 1; // CAL online left 0 - potential isolist pb, keep
    if (RunNumber == 151873) good_run = 1; // same
    if (RunNumber == 151906) good_run = 1; // same
    if (RunNumber == 151907) good_run = 1; // same
    if (RunNumber == 152967) good_run = 0; // CMX offline set to 1 by mistake, remove
    if (RunNumber == 153946) good_run = 0; // SMX HV problems, although SMX online=1, remove
    if (RunNumber == 153985) good_run = 0; // same
//-----------------------------------------------------------------------------
//  3 exceptions from the 2nd part - ask Eric!
//-----------------------------------------------------------------------------
    if (RunNumber == 155324) good_run = 1; // Michael has not looked at it yet, keep
    if (RunNumber == 156370) good_run = 0; // special readout list trig tests, remove
    if (RunNumber == 156371) good_run = 0; // same as above
//-----------------------------------------------------------------------------
//  pre Jan'2003 shutdown runs taken with the PHYSICS_1_04 trigger table
//  (processed around Jan 20 2003, so we can start use them again)
//-----------------------------------------------------------------------------
//      if (RunNumber == 156372) good_run = 0;
//      if (RunNumber == 156401) good_run = 0;
//      if (RunNumber == 156452) good_run = 0;
//      if (RunNumber == 156457) good_run = 0;
//      if (RunNumber == 156460) good_run = 0;
//      if (RunNumber == 156464) good_run = 0;
//-----------------------------------------------------------------------------
//  Jan 13 2003: here the winter'2003 shutdown starts (last physics run 156487)
//-----------------------------------------------------------------------------
    if (RunNumber == 164844) good_run = 0;             // CSL problems
    if (RunNumber == 164870) good_run = 0;             // CSL problems
    if (RunNumber == 164871) good_run = 0;             // CSL problems
    if (RunNumber == 164872) good_run = 0;             // CSL problems
//-----------------------------------------------------------------------------
//  runs "promoted" into GOOD by changing the bits in the database - 
//  make them BAD
//-----------------------------------------------------------------------------
    if (RunNumber == 142558) good_run = 0; // 
    if (RunNumber == 142654) good_run = 0; // 
    if (RunNumber == 144304) good_run = 0; // 
    if (RunNumber == 144624) good_run = 0; // 
    if (RunNumber == 144627) good_run = 0; // 
    if (RunNumber == 144628) good_run = 0; // 
    if (RunNumber == 144713) good_run = 0; // 
    if (RunNumber == 144714) good_run = 0; // 
    if (RunNumber == 145418) good_run = 0; // 

    if ((RunNumber >= 146805) && (RunNumber <= 148157)) good_run = 0; // 

    if (RunNumber == 152680) good_run = 0; // 
    if (RunNumber == 152809) good_run = 0; // 
    if (RunNumber == 153344) good_run = 0; // 

    if (RunNumber == 156372) good_run = 0; // 
    if (RunNumber == 156401) good_run = 0; // 
    if (RunNumber == 156452) good_run = 0; // 
    if (RunNumber == 156457) good_run = 0; // 
    if (RunNumber == 156460) good_run = 0; // 
    if (RunNumber == 156464) good_run = 0; // 
//-----------------------------------------------------------------------------
//  Sep 07 2003: fall'2003 shutdown starts (last physics run: 168892)
//  finally check the trigger table name - make sure it is PHYSICS but not TEST
//-----------------------------------------------------------------------------
  CHECK_TRIGGER_TABLE:
    tt = &rs->TriggerTableName();
    if ((tt->Index("PHYSICS") < 0) || (tt->Index("TEST") >= 0)) {
      good_run = 0;
    }
  }
  return good_run;
}


//_____________________________________________________________________________
Int_t TStnGoodRunList::GoodRunList_DQM_V6(int RunNumber, int RunSection, int Mask) {
  // frozen, runs 141544-183500, DQM V5 good run list
  // total integrated luminosity = ??

  TStnGoodRunList* grl = NULL;

  TString*        tt;
  TStnRunSummary* rs (0);

  int good_run = 0;
					// initialization
  grl = fgGoodRunList;
  
  if ((RunNumber < 141544)) return 0 ;

  rs = grl->GetRunSummary(RunNumber);

  if (rs != 0) {
//-----------------------------------------------------------------------------
// proceed with selecting the good runs: Mario's logic
//-----------------------------------------------------------------------------
    good_run = 1;

    good_run *= rs->ClcStatusBit();
    good_run *= rs->L1tStatusBit();
    good_run *= rs->L2tStatusBit();
    good_run *= rs->L3tStatusBit();
    good_run *= rs->CalStatusBit();
    good_run *= rs->CotStatusBit();
    good_run *= rs->CmuStatusBit();
    good_run *= rs->CmpStatusBit();
    good_run *= rs->SmxStatusBit();

    if (RunNumber < 183500) {
//-----------------------------------------------------------------------------
// 2004-07-07: this is where Mario could get so far
//-----------------------------------------------------------------------------
      good_run *= rs->CalOfflineBit();
      good_run *= rs->CotOfflineBit();
      good_run *= rs->CmuOfflineBit();
      good_run *= rs->CmpOfflineBit();

      if (RunNumber >=150145) {
	good_run *= rs->CmxStatusBit();
	good_run *= rs->CmxOfflineBit();
      }
    }
//-----------------------------------------------------------------------------
// exceptions
//-----------------------------------------------------------------------------
    if (RunNumber == 163438) good_run = 0;
//-----------------------------------------------------------------------------
//  Sep 07 2003: fall'2003 shutdown starts (last physics run: 168892)
//  finally check the trigger table name - make sure it is PHYSICS but not TEST
//-----------------------------------------------------------------------------
    tt = &rs->TriggerTableName();
    if ((tt->Index("PHYSICS") < 0) || (tt->Index("TEST") >= 0)) {
      good_run = 0;
    }
  }
  return good_run;
}

//_____________________________________________________________________________
Int_t TStnGoodRunList::GoodRunList_DQM_V7(int RunNumber, int RunSection, int Mask) {
  // frozen, runs 141544-186598, DQM V7 good run list (up to fall'2004)
  // total integrated luminosity = ??

  static TStnGoodRunList* grl = NULL;

  int electron_flag, muon_flag, silicon_flag, good_run(0);

  TString*        tt;
  TStnRunSummary* rs (0);
					// initialization
  grl = fgGoodRunList;
  
  if      (RunNumber < 141544) {
//-----------------------------------------------------------------------------
//  runs before Feb 04 2002: pre-PHYSICS era
//-----------------------------------------------------------------------------
    return 0 ;
  }
  else if ((RunNumber >= 179096) && (RunNumber <= 182842)) {
//-----------------------------------------------------------------------------
//  runs 179096-182842: running with "compromised" COT 
//-----------------------------------------------------------------------------
    return 0;
  }
  else if (RunNumber > 186598) {
//-----------------------------------------------------------------------------
//  runs 186598+ : future
//-----------------------------------------------------------------------------
    return 0;
  }

  rs = grl->GetRunSummary(RunNumber);

  if (rs != 0) {
//-----------------------------------------------------------------------------
// proceed with selecting the good runs: Mario's logic
//-----------------------------------------------------------------------------
    good_run = 1;
					// CLC
    good_run *= rs->ClcStatusBit();
					// trigger
    good_run *= rs->L1tStatusBit();
    good_run *= rs->L2tStatusBit();
    good_run *= rs->L3tStatusBit();
					// COT (both - online and offline)
    good_run *= rs->CotStatusBit();
    good_run *= rs->CotOfflineBit();
					// calorimeter (online and offline)
    good_run *= rs->CalStatusBit ();
    good_run *= rs->CalOfflineBit();

    electron_flag = grl->ElectronFlag();
    muon_flag     = grl->MuonFlag    ();
    silicon_flag  = grl->SiliconFlag ();
//-----------------------------------------------------------------------------
//  electron bits: 1=CES-only
//                 2=PES-only
//                 3=PES and CES
//-----------------------------------------------------------------------------
    if (electron_flag == 1) {
      good_run *= rs->CesStatusBit ();
      good_run *= rs->CesOfflineBit();
    }
    else if (electron_flag == 2) {
      good_run *= rs->PesStatusBit ();
      good_run *= rs->PesOfflineBit();
    }
    else if (electron_flag == 3) {
      good_run *= rs->CesStatusBit ();
      good_run *= rs->CesOfflineBit();
      good_run *= rs->PesStatusBit ();
      good_run *= rs->PesOfflineBit();
    }
//-----------------------------------------------------------------------------
//  2: CMX: online and offline bits
//-----------------------------------------------------------------------------
    if ((muon_flag == 1) || (muon_flag == 3)) {
      if (RunNumber >= 150145) {
	good_run *= rs->CmxStatusBit ();
	good_run *= rs->CmxOfflineBit();
      }
      else {
	return 0;
      }
    }
//-----------------------------------------------------------------------------
//  3: CMUP: online and offline bits
//-----------------------------------------------------------------------------
    if ((muon_flag == 1) || (muon_flag == 2)) {
      good_run *= rs->CmuStatusBit ();
      good_run *= rs->CmpStatusBit ();
      good_run *= rs->CmuOfflineBit();
      good_run *= rs->CmpOfflineBit();
    }
//-----------------------------------------------------------------------------
//  silicon: online and offline bits
//-----------------------------------------------------------------------------
    if (silicon_flag == 1) {
      good_run *= rs->SvxStatusBit();
      good_run *= rs->IslStatusBit();
      //      good_run *= rs->L00StatusBit();

      good_run *= rs->SvxOfflineBit();
      good_run *= rs->IslOfflineBit();
      //      good_run *= rs->L00OfflineBit();

      if (good_run == 0) return good_run;
    }
//-----------------------------------------------------------------------------
// exceptions
//-----------------------------------------------------------------------------
    if (RunNumber == 163438) good_run = 0; // not sure why
//-----------------------------------------------------------------------------
//  Sep 07 2003: fall'2003 shutdown starts (last physics run: 168892)
//  finally check the trigger table name - make sure it is PHYSICS but not TEST
//-----------------------------------------------------------------------------
    tt = &rs->TriggerTableName();
    if ((tt->Index("PHYSICS") < 0) || (tt->Index("TEST") >= 0)) {
      good_run = 0;
    }
  }
  return good_run;
}


//_____________________________________________________________________________
Int_t TStnGoodRunList::GoodRunList_DQM_V13(int RunNumber,int RunSection, int Mask) {
  // frozen, runs 141544-212133, DQM V13 good run list (up to Feb'2006 shutdown)
  // total integrated luminosity = ??
  // missing so far: partially recovered runs

  static TStnGoodRunList* grl = NULL;

  int electron_flag, muon_flag, silicon_flag, good_run(0);

  TString*        tt;
  TStnRunSummary* rs (0);
					// initialization
  grl = fgGoodRunList;
  
  if      (RunNumber < 141544) {
//-----------------------------------------------------------------------------
//  runs before Feb 04 2002: pre-PHYSICS era
//-----------------------------------------------------------------------------
    return 0 ;
  }
  else if ((RunNumber >= 179096) && (RunNumber <= 182842)) {
//-----------------------------------------------------------------------------
//  runs 179096-182842: running with "compromised" COT 
//-----------------------------------------------------------------------------
    return 0;
  }
  else if (RunNumber > 212133) {
//-----------------------------------------------------------------------------
//  runs 212133+ : future
//-----------------------------------------------------------------------------
    return 0;
  }

  rs = grl->GetRunSummary(RunNumber);

  if (rs != 0) {
//-----------------------------------------------------------------------------
// proceed with selecting the good runs: Mario's logic
//-----------------------------------------------------------------------------
    good_run = 1;
					// CLC
    good_run *= rs->ClcStatusBit();
					// trigger
    good_run *= rs->L1tStatusBit();
    good_run *= rs->L2tStatusBit();
    good_run *= rs->L3tStatusBit();
					// COT (both - online and offline)
    good_run *= rs->CotStatusBit();
    good_run *= rs->CotOfflineBit();
					// calorimeter (online and offline)
    good_run *= rs->CalStatusBit ();
    good_run *= rs->CalOfflineBit();

    electron_flag = grl->ElectronFlag();
    muon_flag     = grl->MuonFlag    ();
    silicon_flag  = grl->SiliconFlag ();
//-----------------------------------------------------------------------------
//  electron bits: 1=CES-only
//                 2=PES-only
//                 3=PES and CES
//-----------------------------------------------------------------------------
    if (electron_flag == 1) {
      good_run *= rs->CesStatusBit ();
      good_run *= rs->CesOfflineBit();
    }
    else if (electron_flag == 2) {
      good_run *= rs->PesStatusBit ();
      good_run *= rs->PesOfflineBit();
    }
    else if (electron_flag == 3) {
      good_run *= rs->CesStatusBit ();
      good_run *= rs->CesOfflineBit();
      good_run *= rs->PesStatusBit ();
      good_run *= rs->PesOfflineBit();
    }
//-----------------------------------------------------------------------------
//  2: CMX: online and offline bits
//     several runs have CMX online bit wrongly set to 'BAD'
//-----------------------------------------------------------------------------
    if ((muon_flag == 1) || (muon_flag == 2)) {
      if (RunNumber >= 150145) {

	if((RunNumber == 183914) || (RunNumber == 183915) || (RunNumber == 183966) ||
	   (RunNumber == 183968) || (RunNumber == 183970) || (RunNumber == 184019) ||
	   (RunNumber == 184020) || (RunNumber == 184021) || (RunNumber == 184028) ||
	   (RunNumber == 184084) || (RunNumber == 184204) || (RunNumber == 184205) ||
	   (RunNumber == 184206) || (RunNumber == 184234) || (RunNumber == 184237) ||
	   (RunNumber == 184240) || (RunNumber == 184243) || (RunNumber == 184289) ||
	   (RunNumber == 184290) || (RunNumber == 184291) || (RunNumber == 184310) ||
	   (RunNumber == 184311) || (RunNumber == 184314) || (RunNumber == 184466) ||
	   (RunNumber == 184467) || (RunNumber == 184469) || (RunNumber == 184516) ||
	   (RunNumber == 184518) || (RunNumber == 185522) || (RunNumber == 185524) ||
	   (RunNumber == 185723) || (RunNumber == 185725) || (RunNumber == 185726) ||
	   (RunNumber == 185969) || (RunNumber == 185970) || (RunNumber == 186040) ||
	   (RunNumber == 186047) || (RunNumber == 186083)) {
	}
	else {
	  // non-exception
	  good_run *= rs->CmxStatusBit ();
	}
	good_run *= rs->CmxOfflineBit();
      }
      else {
	return 0;
      }
    }
//-----------------------------------------------------------------------------
//  3: CMUP: online and offline bits
//-----------------------------------------------------------------------------
    if ((muon_flag == 1) || (muon_flag == 3)) {
      good_run *= rs->CmuStatusBit ();
      good_run *= rs->CmpStatusBit ();
      good_run *= rs->CmuOfflineBit();
      good_run *= rs->CmpOfflineBit();
    }
//-----------------------------------------------------------------------------
//  silicon: online and offline bits
//-----------------------------------------------------------------------------
    if (silicon_flag == 1) {
      good_run *= rs->SvxStatusBit();
      good_run *= rs->IslStatusBit();
      //      good_run *= rs->L00StatusBit();

      good_run *= rs->SvxOfflineBit();
      good_run *= rs->IslOfflineBit();
      //      good_run *= rs->L00OfflineBit();

      if (good_run == 0) return good_run;
    }
//-----------------------------------------------------------------------------
// exceptions
//-----------------------------------------------------------------------------
    if (RunNumber == 163438) good_run = 0; // not sure why
//-----------------------------------------------------------------------------
//  Sep 07 2003: fall'2003 shutdown starts (last physics run: 168892)
//  finally check the trigger table name - make sure it is PHYSICS but not TEST
//  Skip diffractive runs taken at low luminosity:
//               209865, 209866, 211058, 211073, 211079
//-----------------------------------------------------------------------------
    tt = &rs->TriggerTableName();
    if ((tt->Index("PHYSICS" ) <  0) || 
	(tt->Index("TEST"    ) >= 0) || 
	(tt->Index("DIFFLOWL") >= 0)    ) {
      good_run = 0;
    }
  }
  return good_run;
}


//_____________________________________________________________________________
Int_t TStnGoodRunList::GoodRunList_DQM_V27(int RunNumber,int RunSection, int Mask) {
  // frozen, runs 141544-272114, DQM V27 good run list (up to spring 2009)
  // total integrated luminosity = ??

  TStnGoodRunList* grl(0);

  int electron_flag, muon_flag, silicon_flag, good_run(0);

  TString*        tt;
  TStnRunSummary* rs (0);
					// initialization
  grl = fgGoodRunList;
//-----------------------------------------------------------------------------
//  runs before Feb 04 2002: pre-PHYSICS era
//-----------------------------------------------------------------------------
  if (RunNumber <  141544)                                  return 0;
//-----------------------------------------------------------------------------
// not sure why run 161438 is excluded, but it is...
//-----------------------------------------------------------------------------
  if (RunNumber == 163438)                                  return 0;
//-----------------------------------------------------------------------------
// runs 179096-182842: "compromised" COT 
//-----------------------------------------------------------------------------
  if ((RunNumber >= 179096) && (RunNumber <= 182842))       return 0;
//------------------------------------------------------------------------------
// COT SL2 problems
//-----------------------------------------------------------------------------
  if ( (RunNumber == 264810) || 
       (RunNumber == 264811) || 
      ((RunNumber >= 264854) && (RunNumber <= 264991)) )    return 0;
//-----------------------------------------------------------------------------
//  runs 271047+ : future
//-----------------------------------------------------------------------------
  if (RunNumber > 271047)                                   return 0;

  rs = grl->GetRunSummary(RunNumber);
  if (rs == 0)                                              return 0;
//-----------------------------------------------------------------------------
// in case of partially recovered runs can't believe to the bit settings
// partially recovered runs in my are supposed to be ALL GOOD for the good 
// run section ranges 
// start from checking partially recovered runs: for them we can't rely on
// the bit settings, if NRsRanges > 0, the run is partially recovered
// consider partially recovered run as good within the range of good run 
// sections
// these are physics runs by construction - a recovery has been attempted
// but the run section has to be in one of the good ranges
//-----------------------------------------------------------------------------
  int rr  = 0;
  int nrr = rs->NRsRanges();
  good_run = 0;
  for (int irr=0; irr<nrr; irr++) {
    if (RunSection < 0) {
				// choose the first runsection range....
      rr       = irr;
      good_run = 1;
      break;
    }
    else if ((RunSection >= rs->LowRS(irr)) && (RunSection <= rs->HighRS(irr))) {
      rr       = irr;
      good_run = 1;
      break;
    }
  }
  if (good_run == 0) return 0;
//-----------------------------------------------------------------------------
// now proceed with selecting the good runs: Mario's logic
//-----------------------------------------------------------------------------
  good_run = 1;
					// CLC
  good_run *= rs->ClcStatusBit(rr);
					// trigger
  good_run *= rs->L1tStatusBit(rr);
  good_run *= rs->L2tStatusBit(rr);
  good_run *= rs->L3tStatusBit(rr);
					// COT (both - online and offline)
  good_run *= rs->CotStatusBit (rr);
  good_run *= rs->CotOfflineBit(rr);
					// calorimeter (online and offline)
  good_run *= rs->CalStatusBit (rr);
  good_run *= rs->CalOfflineBit(rr);

  electron_flag = grl->ElectronFlag();
  muon_flag     = grl->MuonFlag    ();
  silicon_flag  = grl->SiliconFlag ();
//-----------------------------------------------------------------------------
//  electron bits: 1=CES-only
//                 2=PES-only
//                 3=PES and CES
//-----------------------------------------------------------------------------
  if (electron_flag == 1) {
    good_run *= rs->CesStatusBit (rr);
    good_run *= rs->CesOfflineBit(rr);
  }
  else if (electron_flag == 2) {
    good_run *= rs->PesStatusBit (rr);
    good_run *= rs->PesOfflineBit(rr);
  }
  else if (electron_flag == 3) {
    good_run *= rs->CesStatusBit (rr);
    good_run *= rs->CesOfflineBit(rr);
    good_run *= rs->PesStatusBit (rr);
    good_run *= rs->PesOfflineBit(rr);
  }
//-----------------------------------------------------------------------------
// 2: CMX: online and offline bits. CMX operational starting from run 150145
//    muon_flag=1: CMUP+CMX, 3: CMX-only 
//-----------------------------------------------------------------------------
  if ((muon_flag == 1) || (muon_flag == 2)) {
    if (RunNumber < 150145)                                 return 0;
//-----------------------------------------------------------------------------
// runs for which CMX was marked as bad ONLINE, but it should be ignored
// for these runs check only offline bit
//-----------------------------------------------------------------------------
    if((RunNumber == 183914) || (RunNumber == 183915) || (RunNumber == 183966) ||
       (RunNumber == 183968) || (RunNumber == 183970) || (RunNumber == 184019) ||
       (RunNumber == 184020) || (RunNumber == 184021) || (RunNumber == 184028) ||
       (RunNumber == 184084) || (RunNumber == 184204) || (RunNumber == 184205) ||
       (RunNumber == 184206) || (RunNumber == 184234) || (RunNumber == 184237) ||
       (RunNumber == 184240) || (RunNumber == 184243) || (RunNumber == 184289) ||
       (RunNumber == 184290) || (RunNumber == 184291) || (RunNumber == 184310) ||
       (RunNumber == 184311) || (RunNumber == 184314) || (RunNumber == 184466) ||
       (RunNumber == 184467) || (RunNumber == 184469) || (RunNumber == 184516) ||
       (RunNumber == 184518) || (RunNumber == 185522) || (RunNumber == 185524) ||
       (RunNumber == 185723) || (RunNumber == 185725) || (RunNumber == 185726) ||
       (RunNumber == 185969) || (RunNumber == 185970) || (RunNumber == 186040) ||
       (RunNumber == 186047) || (RunNumber == 186083)) {
      good_run *= rs->CmxOfflineBit(rr);
    }
    else {
      good_run *= rs->CmxStatusBit (rr);
      good_run *= rs->CmxOfflineBit(rr);
    }
  }
//-----------------------------------------------------------------------------
//  3: CMUP: online and offline bits
//-----------------------------------------------------------------------------
  if ((muon_flag == 1) || (muon_flag == 2)) {
    good_run *= rs->CmuStatusBit (rr);
    good_run *= rs->CmpStatusBit (rr);
    good_run *= rs->CmuOfflineBit(rr);
    good_run *= rs->CmpOfflineBit(rr);
  }
//-----------------------------------------------------------------------------
//  silicon: online and offline bits
//-----------------------------------------------------------------------------
  if (silicon_flag == 1) {
//-----------------------------------------------------------------------------
// SVX low efficiency problem
//-----------------------------------------------------------------------------
    if ((RunNumber >= 262865) && (RunNumber <= 263615)) {
      good_run = 0;
    }
    else {
      good_run *= rs->SvxStatusBit(rr);
      good_run *= rs->IslStatusBit(rr);
      //      good_run *= rs->L00StatusBit();

      good_run *= rs->SvxOfflineBit(rr);
      good_run *= rs->IslOfflineBit(rr);
      //      good_run *= rs->L00OfflineBit();
    }
  }

  if (good_run != 0) {
//-----------------------------------------------------------------------------
//  finally check the trigger table name - make sure it is PHYSICS but not TEST
//-----------------------------------------------------------------------------
    tt = &rs->TriggerTableName();
    if ((tt->Index("PHYSICS") < 0) || (tt->Index("TEST") >= 0)) {
      good_run = 0;
    }
  }

  return good_run;
}


//_____________________________________________________________________________
Int_t TStnGoodRunList::GoodRunList_DQM_V32(int RunNumber,int RunSection, int Mask) {
  // frozen, runs 141544-284843, DQM V32 good run list (up to Oct 25 2009)
  // total integrated luminosity = ??

  TStnGoodRunList* grl(0);

  int electron_flag, muon_flag, silicon_flag, good_run(0);

  TString*        tt;
  TStnRunSummary* rs (0);
					// initialization
  grl = fgGoodRunList;
//-----------------------------------------------------------------------------
//  runs before Feb 04 2002: pre-PHYSICS era
//-----------------------------------------------------------------------------
  if (RunNumber <  141544)                                  return 0;
//-----------------------------------------------------------------------------
// not sure why run 161438 is excluded, but it is...
//-----------------------------------------------------------------------------
  if (RunNumber == 163438)                                  return 0;
//-----------------------------------------------------------------------------
// runs 179096-182842: "compromised" COT 
//-----------------------------------------------------------------------------
  if ((RunNumber >= 179096) && (RunNumber <= 182842))       return 0;
//------------------------------------------------------------------------------
// COT SL2 problems
//-----------------------------------------------------------------------------
  if ( (RunNumber == 264810) || 
       (RunNumber == 264811) || 
      ((RunNumber >= 264854) && (RunNumber <= 264991)) )    return 0;
//-----------------------------------------------------------------------------
// run 267255: L1 bad
//-----------------------------------------------------------------------------
  if (RunNumber == 267255)                                  return 0;
//-----------------------------------------------------------------------------
// run 272470: no minbias trigger
// run 273251: bad CMU, they also say it has not been processed
//-----------------------------------------------------------------------------
  if ((RunNumber == 272470) || (RunNumber == 273251))       return 0;
//-----------------------------------------------------------------------------
//  runs 284843+ : future
//-----------------------------------------------------------------------------
  if (RunNumber > 284843)                                   return 0;

  rs = grl->GetRunSummary(RunNumber);
  if (rs == 0)                                              return 0;
//-----------------------------------------------------------------------------
// in case of partially recovered runs can't believe to the bit settings
// partially recovered runs in my are supposed to be ALL GOOD for the good 
// run section ranges 
// start from checking partially recovered runs: for them we can't rely on
// the bit settings, if NRsRanges > 0, the run is partially recovered
// consider partially recovered run as good within the range of good run 
// sections
// these are physics runs by construction - a recovery has been attempted
// but the run section has to be in one of the good ranges
//-----------------------------------------------------------------------------
  int rr  = 0;
  int nrr = rs->NRsRanges();
  good_run = 0;
  for (int irr=0; irr<nrr; irr++) {
    if (RunSection < 0) {
				// choose the first runsection range....
      rr       = irr;
      good_run = 1;
      break;
    }
    else if ((RunSection >= rs->LowRS(irr)) && (RunSection <= rs->HighRS(irr))) {
      rr       = irr;
      good_run = 1;
      break;
    }
  }
  if (good_run == 0) return 0;
//-----------------------------------------------------------------------------
// now proceed with selecting the good runs: Mario's logic
//-----------------------------------------------------------------------------
  good_run = 1;
					// CLC
  good_run *= rs->ClcStatusBit(rr);
					// trigger
  good_run *= rs->L1tStatusBit(rr);
  good_run *= rs->L2tStatusBit(rr);
  good_run *= rs->L3tStatusBit(rr);
					// COT (both - online and offline)
  good_run *= rs->CotStatusBit (rr);
  good_run *= rs->CotOfflineBit(rr);
					// calorimeter (online and offline)
  good_run *= rs->CalStatusBit (rr);
  good_run *= rs->CalOfflineBit(rr);

  electron_flag = grl->ElectronFlag();
  muon_flag     = grl->MuonFlag    ();
  silicon_flag  = grl->SiliconFlag ();
//-----------------------------------------------------------------------------
//  electron bits: 1=CES-only
//                 2=PES-only
//                 3=PES and CES
//-----------------------------------------------------------------------------
  if (electron_flag == 1) {
    good_run *= rs->SmxStatusBit (rr);
    //    good_run *= rs->CesOfflineBit(rr);
  }
//   else if (electron_flag == 2) {
//     good_run *= rs->PesStatusBit (rr);
//     good_run *= rs->PesOfflineBit(rr);
//   }
  else if (electron_flag == 3) {
    good_run *= rs->SmxStatusBit (rr);
//     good_run *= rs->CesOfflineBit(rr);
//     good_run *= rs->PesStatusBit (rr);
//     good_run *= rs->PesOfflineBit(rr);
  }
//-----------------------------------------------------------------------------
// 2: CMX: online and offline bits. CMX operational starting from run 150145
//    muon_flag=1: CMUP+CMX, 3: CMX-only 
//-----------------------------------------------------------------------------
  if ((muon_flag == 1) || (muon_flag == 3)) {
    if (RunNumber < 150145)                                 return 0;
//-----------------------------------------------------------------------------
// runs for which CMX was marked as bad ONLINE, but it should be ignored
// for these runs check only offline bit
//-----------------------------------------------------------------------------
    if((RunNumber == 183914) || (RunNumber == 183915) || (RunNumber == 183966) ||
       (RunNumber == 183968) || (RunNumber == 183970) || (RunNumber == 184019) ||
       (RunNumber == 184020) || (RunNumber == 184021) || (RunNumber == 184028) ||
       (RunNumber == 184084) || (RunNumber == 184204) || (RunNumber == 184205) ||
       (RunNumber == 184206) || (RunNumber == 184234) || (RunNumber == 184237) ||
       (RunNumber == 184240) || (RunNumber == 184243) || (RunNumber == 184289) ||
       (RunNumber == 184290) || (RunNumber == 184291) || (RunNumber == 184310) ||
       (RunNumber == 184311) || (RunNumber == 184314) || (RunNumber == 184466) ||
       (RunNumber == 184467) || (RunNumber == 184469) || (RunNumber == 184516) ||
       (RunNumber == 184518) || (RunNumber == 185522) || (RunNumber == 185524) ||
       (RunNumber == 185723) || (RunNumber == 185725) || (RunNumber == 185726) ||
       (RunNumber == 185969) || (RunNumber == 185970) || (RunNumber == 186040) ||
       (RunNumber == 186047) || (RunNumber == 186083)) {
      good_run *= rs->CmxOfflineBit(rr);
    }
    else {
      good_run *= rs->CmxStatusBit (rr);
      good_run *= rs->CmxOfflineBit(rr);
    }
  }
//-----------------------------------------------------------------------------
//  3: CMUP: online and offline bits +IMU
//-----------------------------------------------------------------------------
  if ((muon_flag == 1) || (muon_flag == 2)) {
    good_run *= rs->CmuStatusBit (rr);
    good_run *= rs->CmpStatusBit (rr);
    good_run *= rs->ImuStatusBit (rr);

    good_run *= rs->CmuOfflineBit(rr);
    good_run *= rs->CmpOfflineBit(rr);
    good_run *= rs->ImuOfflineBit(rr);
  }
//-----------------------------------------------------------------------------
//  silicon: online and offline bits
//-----------------------------------------------------------------------------
  if (silicon_flag == 1) {
//-----------------------------------------------------------------------------
// SVX low efficiency problem
//-----------------------------------------------------------------------------
    if      ((RunNumber >= 262865) && (RunNumber <= 263615)) {
      good_run = 0;
    }
    else if ((RunNumber == 267335) || 
	     (RunNumber == 268978) || 
	     (RunNumber == 277340)    ) {
      good_run = 0;
    }
    else {
      good_run *= rs->SvxStatusBit(rr);
      good_run *= rs->SvxOfflineBit(rr);

      if ((RunNumber >= 241665) && (RunNumber <= 246231)) {
				        // ignore ISL bit for period 13 (241665-246231)
      }
      else {
	good_run *= rs->IslStatusBit (rr);
	good_run *= rs->IslOfflineBit(rr);
      }
      //      good_run *= rs->L00StatusBit();
      //      good_run *= rs->L00OfflineBit();
    }
  }

  if (good_run != 0) {
//-----------------------------------------------------------------------------
//  finally check the trigger table name - make sure it is PHYSICS but not TEST
//-----------------------------------------------------------------------------
    tt = &rs->TriggerTableName();
    if ((tt->Index("PHYSICS") < 0) || (tt->Index("TEST") >= 0)) {
      good_run = 0;
    }
  }

  return good_run;
}


//_____________________________________________________________________________
Int_t TStnGoodRunList::GoodRunList_DQM_V34(int RunNumber, int RunSection, int Mask) {
  // frozen, runs 141544-289197, DQM V34 good run list (up to Feb 25 2010)
  // total integrated luminosity = ??
  // if Mask != 0 then it defines electron, muon and silicon bits

  TStnGoodRunList* grl(0);

  int electron_flag, muon_flag, silicon_flag, good_run(0);

  TString*        tt;
  TStnRunSummary* rs (0);
					// initialization
  grl = fgGoodRunList;
//-----------------------------------------------------------------------------
//  runs before Feb 04 2002: pre-PHYSICS era
//-----------------------------------------------------------------------------
  if (RunNumber <  141544)                                  return 0;
//-----------------------------------------------------------------------------
// not sure why run 161438 is excluded, but it is...
//-----------------------------------------------------------------------------
  if (RunNumber == 163438)                                  return 0;
//-----------------------------------------------------------------------------
// runs 179096-182842: "compromised" COT 
//-----------------------------------------------------------------------------
  if ((RunNumber >= 179096) && (RunNumber <= 182842))       return 0;
//------------------------------------------------------------------------------
// COT SL2 problems
//-----------------------------------------------------------------------------
  if ( (RunNumber == 264810) || 
       (RunNumber == 264811) || 
      ((RunNumber >= 264854) && (RunNumber <= 264991)) )    return 0;
//-----------------------------------------------------------------------------
// run 267255: L1 bad
//-----------------------------------------------------------------------------
  if (RunNumber == 267255)                                  return 0;
//-----------------------------------------------------------------------------
// run 272470: no minbias trigger
// run 273251: bad CMU, they also say it has not been processed
//-----------------------------------------------------------------------------
  if ((RunNumber == 272470) || (RunNumber == 273251))       return 0;
//-----------------------------------------------------------------------------
//  runs 289197+ : future
//-----------------------------------------------------------------------------
  if (RunNumber > 2889197)                                  return 0;

  rs = grl->GetRunSummary(RunNumber);
  if (rs == 0)                                              return 0;
//-----------------------------------------------------------------------------
// in case of partially recovered runs can't believe to the bit settings
// partially recovered runs in my are supposed to be ALL GOOD for the good 
// run section ranges 
// start from checking partially recovered runs: for them we can't rely on
// the bit settings, if NRsRanges > 0, the run is partially recovered
// consider partially recovered run as good within the range of good run 
// sections
// these are physics runs by construction - a recovery has been attempted
// but the run section has to be in one of the good ranges
//-----------------------------------------------------------------------------
  int rr  = 0;
  int nrr = rs->NRsRanges();
  good_run = 0;
  for (int irr=0; irr<nrr; irr++) {
    if (RunSection < 0) {
				// choose the first runsection range....
      rr       = irr;
      good_run = 1;
      break;
    }
    else if ((RunSection >= rs->LowRS(irr)) && (RunSection <= rs->HighRS(irr))) {
      rr       = irr;
      good_run = 1;
      break;
    }
  }
  if (good_run == 0) return 0;
//-----------------------------------------------------------------------------
// now proceed with selecting the good runs: Mario's logic
//-----------------------------------------------------------------------------
  good_run = 1;
					// CLC
  good_run *= rs->ClcStatusBit(rr);
					// trigger
  good_run *= rs->L1tStatusBit(rr);
  good_run *= rs->L2tStatusBit(rr);
  good_run *= rs->L3tStatusBit(rr);
					// COT (both - online and offline)
  good_run *= rs->CotStatusBit (rr);
  good_run *= rs->CotOfflineBit(rr);
					// calorimeter (online and offline)
  good_run *= rs->CalStatusBit (rr);
  good_run *= rs->CalOfflineBit(rr);

  if (Mask == 0) {
    electron_flag = grl->ElectronFlag();
    muon_flag     = grl->MuonFlag    ();
    silicon_flag  = grl->SiliconFlag ();
  }
  else {
    electron_flag = Mask & 0x1;
    muon_flag     = Mask & 0x2;
    silicon_flag  = Mask & 0x4;
  }
//-----------------------------------------------------------------------------
//  electron bits: 1=CES-only
//                 2=PES-only
//                 3=PES and CES
//-----------------------------------------------------------------------------
  if (electron_flag == 1) {
    good_run *= rs->SmxStatusBit (rr);
    //    good_run *= rs->CesOfflineBit(rr);
  }
//   else if (electron_flag == 2) {
//     good_run *= rs->PesStatusBit (rr);
//     good_run *= rs->PesOfflineBit(rr);
//   }
  else if (electron_flag == 3) {
    good_run *= rs->SmxStatusBit (rr);
//     good_run *= rs->CesOfflineBit(rr);
//     good_run *= rs->PesStatusBit (rr);
//     good_run *= rs->PesOfflineBit(rr);
  }
//-----------------------------------------------------------------------------
// 2: CMX: online and offline bits. CMX operational starting from run 150145
//    muon_flag=1: CMUP+CMX, 3: CMX-only 
//-----------------------------------------------------------------------------
  if ((muon_flag == 1) || (muon_flag == 3)) {
    if (RunNumber < 150145)                                 return 0;
//-----------------------------------------------------------------------------
// runs for which CMX was marked as bad ONLINE, but it should be ignored
// for these runs check only offline bit
//-----------------------------------------------------------------------------
    if((RunNumber == 183914) || (RunNumber == 183915) || (RunNumber == 183966) ||
       (RunNumber == 183968) || (RunNumber == 183970) || (RunNumber == 184019) ||
       (RunNumber == 184020) || (RunNumber == 184021) || (RunNumber == 184028) ||
       (RunNumber == 184084) || (RunNumber == 184204) || (RunNumber == 184205) ||
       (RunNumber == 184206) || (RunNumber == 184234) || (RunNumber == 184237) ||
       (RunNumber == 184240) || (RunNumber == 184243) || (RunNumber == 184289) ||
       (RunNumber == 184290) || (RunNumber == 184291) || (RunNumber == 184310) ||
       (RunNumber == 184311) || (RunNumber == 184314) || (RunNumber == 184466) ||
       (RunNumber == 184467) || (RunNumber == 184469) || (RunNumber == 184516) ||
       (RunNumber == 184518) || (RunNumber == 185522) || (RunNumber == 185524) ||
       (RunNumber == 185723) || (RunNumber == 185725) || (RunNumber == 185726) ||
       (RunNumber == 185969) || (RunNumber == 185970) || (RunNumber == 186040) ||
       (RunNumber == 186047) || (RunNumber == 186083)) {
      good_run *= rs->CmxOfflineBit(rr);
    }
    else {
      good_run *= rs->CmxStatusBit (rr);
      good_run *= rs->CmxOfflineBit(rr);
    }
  }
//-----------------------------------------------------------------------------
//  3: CMUP: online and offline bits +IMU
//-----------------------------------------------------------------------------
  if ((muon_flag == 1) || (muon_flag == 2)) {
    good_run *= rs->CmuStatusBit (rr);
    good_run *= rs->CmpStatusBit (rr);
    good_run *= rs->ImuStatusBit (rr);

    good_run *= rs->CmuOfflineBit(rr);
    good_run *= rs->CmpOfflineBit(rr);
    good_run *= rs->ImuOfflineBit(rr);
  }
//-----------------------------------------------------------------------------
//  silicon: online and offline bits
//-----------------------------------------------------------------------------
  if (silicon_flag == 1) {
//-----------------------------------------------------------------------------
// SVX low efficiency problem
//-----------------------------------------------------------------------------
    if      ((RunNumber >= 262865) && (RunNumber <= 263615)) {
      good_run = 0;
    }
    else if ((RunNumber == 267335) || 
	     (RunNumber == 268978) || 
	     (RunNumber == 277340)    ) {
      good_run = 0;
    }
    else {
      good_run *= rs->SvxStatusBit(rr);
      good_run *= rs->SvxOfflineBit(rr);

      if ((RunNumber >= 241665) && (RunNumber <= 246231)) {
				        // ignore ISL bit for period 13 (241665-246231)
      }
      else {
	good_run *= rs->IslStatusBit (rr);
	good_run *= rs->IslOfflineBit(rr);
      }
      //      good_run *= rs->L00StatusBit();
      //      good_run *= rs->L00OfflineBit();
    }
  }

  if (good_run != 0) {
//-----------------------------------------------------------------------------
//  finally check the trigger table name - make sure it is PHYSICS but not TEST
//-----------------------------------------------------------------------------
    tt = &rs->TriggerTableName();
    if ((tt->Index("PHYSICS") < 0) || (tt->Index("TEST") >= 0)) {
      good_run = 0;
    }
  }

  return good_run;
}


//_____________________________________________________________________________
Int_t TStnGoodRunList::GoodRunList_MC_1001(int RunNumber, int RunSection, int Mask) {
  // the reason for having this good run list is that 
  // mcProduction/base/runlist_v7goodrunlist_1001 used to generate bulk of 
  // 5.3.3 datasets included set of runs without the silicon, which, however,
  // had silicon beam positions defined and set to (0.1,0.5)
  // total integrated luminosity = ??
  // MC statistics is large enough, always use MC run set which is a subset 
  // of the data run set


  int good_run;

  good_run = GoodRunList_DQM_V34(RunNumber,RunSection,Mask);

  if       (RunNumber  < 141572)                           good_run = 0;
  else if ((RunNumber >= 150951) && (RunNumber <= 151240)) good_run = 0;
  else if ((RunNumber >= 154049) && (RunNumber <= 154070)) good_run = 0;
  else if ((RunNumber >= 154575) && (RunNumber <= 154613)) good_run = 0;
  else if ((RunNumber >= 154710) && (RunNumber <= 154861)) good_run = 0;
  else if  (RunNumber == 155155)                           good_run = 0;
  else if ((RunNumber >= 158644) && (RunNumber <= 158799)) good_run = 0;
  else if  (RunNumber == 160259)                           good_run = 0;
  else if  (RunNumber == 167770)                           good_run = 0;
  else if  (RunNumber == 168125)                           good_run = 0;
  else if ((RunNumber >= 175585) && (RunNumber <= 176646)) good_run = 0;
  else if ((RunNumber >= 177140) && (RunNumber <= 177158)) good_run = 0;
  else if ((RunNumber >= 177620) && (RunNumber <= 177621)) good_run = 0;
  else if  (RunNumber == 178254)                           good_run = 0;
  else if  (RunNumber == 178388)                           good_run = 0;
  else if ((RunNumber >= 179234) && (RunNumber <= 179463)) good_run = 0;
  else if ((RunNumber >= 179578) && (RunNumber <= 179621)) good_run = 0;
  else if ((RunNumber >= 180958) && (RunNumber <= 180965)) good_run = 0;
  else if ((RunNumber >= 181801) && (RunNumber <= 181880)) good_run = 0;
  else if  (RunNumber == 182111)                           good_run = 0;
  else if  (RunNumber == 182134)                           good_run = 0;
  else if  (RunNumber == 182217)                           good_run = 0;
  else if ((RunNumber >= 182335) && (RunNumber <= 182337)) good_run = 0;
  else if ((RunNumber >= 182695) && (RunNumber <= 182719)) good_run = 0;
  else if ((RunNumber >= 182889) && (RunNumber <= 182898)) good_run = 0;
  else if ((RunNumber >= 183510) && (RunNumber <= 183511)) good_run = 0;

  return good_run;
}


//_____________________________________________________________________________
int TStnGoodRunList::DefaultGoodRunRoutine(Int_t RunNumber, int RunSection, int Mask) {
  // to start with assume that the constants for a given run are cached,
  // currently the interface is ahead of the implementation

  int  rc(0);

//    static int              nlines          = 0;
//    static TStnGoodRunList* listOfGoodRuns  = NULL;
//    static TFolder*         folder          = NULL;
//    static int              print_level     = 0;

//    TStnRunSummary* rs;
//    Int_t           rmin, rmax, nr;
//    const char*     text;

//    const    char*  GOOD        = "GOOD";
//    const    char*  CHECKED_BAD = "CHECKED, BAD";
//    const    char*  LOW_LUMI    = "LOW LUMI";
//    const    char*  BAD_GRL     = "GRL: BAD";
//    const    char*  UNCHECKED   = "UNCHECKED";

//  //-----------------------------------------------------------------------------
//  // check for the good Run List, first call: RunNumber = print level
//  //-----------------------------------------------------------------------------
//    if (RunNumber <= 0) {
//      folder         = (TFolder*) gROOT->GetRootFolder()->FindObject("Ana");
//      listOfGoodRuns = (TStnGoodRunList*) folder->FindObject("GoodRunList");
//      print_level    = -RunNumber;
//    }
//  //-----------------------------------------------------------------------------
//  // if list of good runs is not defined, everything is good
//  //-----------------------------------------------------------------------------
//    TStnDBManager* dbm = TStnDBManager::Instance();

//    rs       = (TStnRunSummary*) dbm->GetTable("RunSummary");

//    if (! listOfGoodRuns) {
//      rc   = 1;
//      text = GOOD;
//                                                              goto END;
//    }
//  //-----------------------------------------------------------------------------
//  // make sure run is GOOD and ON TAPE luminosity is above the threshold
//  // "RUN IS GOOD" =
//  // - either GoodRunBit is set
//  // - or run not checked
//  //-----------------------------------------------------------------------------
//    if (rs->IsChecked()) {
//  //-----------------------------------------------------------------------------
//  // checked run has to have GOOD RUN bit set to 1
//  //-----------------------------------------------------------------------------
//      if (rs->GoodRunStatus() == 0) {
//        text = CHECKED_BAD;
//        rc   = 0;
//                                                              goto END;
//      }
//    }
//    else {
//  //-----------------------------------------------------------------------------
//  // not checked run : see if we deal with them
//  //-----------------------------------------------------------------------------
//      if (listOfGoodRuns->UseUncheckedRuns() == 0) {
//        text = UNCHECKED;
//        rc   = 0;
//                                                              goto END;
//      }
//    }
//  //-----------------------------------------------------------------------------
//  // check run luminosity
//  //-----------------------------------------------------------------------------
//    if (rs->LumiTape() < listOfGoodRuns->MinGoodRunLum()) {
//      rc   = 0;
//      text = LOW_LUMI;
//                                                              goto END;
//    }
//  //-----------------------------------------------------------------------------
//  // check if the run has not been marked as bad by the user who defined GRL
//  //-----------------------------------------------------------------------------
//    text  = BAD_GRL;
//    nr = listOfGoodRuns->NRanges();
//    rc = 0;
//    for (int i=0; i<nr; i++) {
//      rmin = listOfGoodRuns->MinRunNumber(i);
//      rmax = listOfGoodRuns->MaxRunNumber(i);
//      if (RunNumber < rmin) {
//  //-----------------------------------------------------------------------------
//  // assume that the run ranges are ordered according to accending run number
//  //-----------------------------------------------------------------------------
//        break;
//      }
//      if ((RunNumber >= rmin) && (RunNumber <= rmax)) {
//        rc    = 1;
//        text  = GOOD;
//        break;
//      }
//    }
//  //-----------------------------------------------------------------------------
//  // make sure nothing gets unnoticed
//  //-----------------------------------------------------------------------------
//   END:
//    if (print_level == 10) {
//      if (nlines == 0) {
//        printf("--------------------------------------------------");
//        printf("----------------------\n");
//        printf("      run   rc status                   ");
//        printf(" L(TeV)     L(live)   L(offline)\n");
//        printf("---------------------------------------------------");
//        printf("---------------------\n");
//      }
//      printf(" %8i  %3i %-20s  %10.3f  %10.3f   %10.3f\n",
//  	   rs->RunNumber(),
//  	   rc,
//  	   text,
//  	   rs->LumiTev(),
//  	   rs->LumiTape(),
//  	   rs->OfflineLumiRS());
//      nlines++;
//      if (nlines == 50) nlines = 0;
//    }
  return rc;
}


//-----------------------------------------------------------------------------
void TStnGoodRunList::Clear(const char* Opt) {
}

//-----------------------------------------------------------------------------
void TStnGoodRunList::Print(const char* Opt) const {
}

