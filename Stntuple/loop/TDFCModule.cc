///////////////////////////////////////////////////////////////////////////////
// this module runs on a single STNTUPLE file and collects information
// necessary to make a Data File Catalog (DFC) entry
// fPrintLevel =    1 :
//                  2 : 1 + run-level catalog of each file
//                  3 : getFileInfo-style printout for grabber
//             >   10 : also print commands to be executed for cataloging
//
// fExecCommands =  1 : execute move command
//               = 10 : execute cataloging command
//               = 11 : execute both commands
///////////////////////////////////////////////////////////////////////////////
#include "iostream"
using namespace std;
#include "TSystem.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TText.h"
#include "TChain.h"

#include "Stntuple/loop/TStnAna.hh"
#include "Stntuple/loop/TStnInputModule.hh"
#include "Stntuple/obj/TStnHeaderBlock.hh"

#include "Stntuple/loop/TDFCModule.hh"

ClassImp(TDFCModule)
//_____________________________________________________________________________
TDFCModule::TDFCModule(const char* name, const char* title):
  TStnModule(name,title),
  fDatasetID("UNKNWN"),
  fBook     ("CDFPBOOK"),
  fDbID     ("production_file_catalog_write"),
  fPrintOpt ("banner")
{
//-----------------------------------------------------------------------------
// by default print run-level catalog of each file (fPrintLevel=2)
//-----------------------------------------------------------------------------
  fMinRunNumber     = 1000000;
  fMinSectionNumber = 1000000;
  fMinEventNumber   = 100000000;
  fMaxRunNumber     = -1;
  fMaxEventNumber   = -1;
  fFile             = 0;
  fFileSize         = 0.;
  fFileDate         = TString("");
  fNewFile          = 0;
  fNEvents          = 0;
  fExecCommands     = 0;
  fPrintLevel       = 2;
  fNFileset         = 0;
  fReturnCode       = 0;
  fOutputDir        = 0;
  fCurrentRunRecord = 0;
  fRunRecord.clear();
}

//_____________________________________________________________________________
TDFCModule::~TDFCModule()
{
  int nr = fRunRecord.size();
  for (int i=0; i<nr; i++) {
    delete fRunRecord[i];
  }

  delete fNewFile;
  delete fOutputDir;
}

//_____________________________________________________________________________
int TDFCModule::BeginJob()
{
  // book histograms, this code accomodates for run number < 200000
  TChain* chain = GetAna()->GetInputModule()->GetChain();
  fFile         = chain->GetFile();
  fNEvents      = 0;

  char cmd[200],fn[200]="";
  sprintf(cmd,"echo %s |  awk -F / ' {print $NF}'",fFile->GetName());
  FILE* f = gSystem->OpenPipe(cmd,"r");
  fscanf(f,"%s",fn);
  fclose(f);

  fFileName = fn;
  fFileSize = fFile->GetSize()/1000000.;
  fFileDate = TString(fFile->GetModificationDate().AsSQLString());

  return 0;
}


//_____________________________________________________________________________
int TDFCModule::BeginRun() {
//---------------------------------------------------------------------------
// Figure out whether we know this run, otherwise install a new record
// Try to find this run in our filing system
//---------------------------------------------------------------------------
  int  iRun  = GetHeaderBlock()->RunNumber();
  bool found = false;
  int nr = fRunRecord.size();
  for (int i=0; i<nr; i++) {
    if (fRunRecord[i]->fNumber == iRun) {
      // Okay this is our run
      fCurrentRunRecord = fRunRecord[i];
      found = true;
    }
  }
  // This run was not yet found.... create new record
  if (! found) {

    fCurrentRunRecord = new RunRecord_t;
    fCurrentRunRecord->fNumber         = -1;
    fCurrentRunRecord->fNEvents        =  0;
    fCurrentRunRecord->fLoEvent        = -1;
    fCurrentRunRecord->fHiEvent        = -1;
    fCurrentRunRecord->fListOfRunSections.clear();

    fRunRecord.push_back(fCurrentRunRecord);
  }

  return 0;
}

//_____________________________________________________________________________
int TDFCModule::Event(int ientry)
{
  TStnHeaderBlock* header = GetHeaderBlock();

  int rn = header->RunNumber    ();
  int ev = header->EventNumber  ();
  int rs = header->SectionNumber();

  if (fCurrentRunRecord->fNumber < 0) {
    fCurrentRunRecord->fNumber  = rn;
    fCurrentRunRecord->fLoEvent = ev;
    fCurrentRunRecord->fHiEvent = ev;
  }
  else if (ev > fCurrentRunRecord->fHiEvent)
    fCurrentRunRecord->fHiEvent = ev;
  else if (ev < fCurrentRunRecord->fLoEvent)
    fCurrentRunRecord->fLoEvent = ev;

  fCurrentRunRecord->fNEvents++;
  
  if (rn < fMinRunNumber) {
    fMinRunNumber     = rn;
    fMinSectionNumber = 1000000;
    fMinEventNumber   = 1000000000;
  }

  if (rn == fMinRunNumber) {
    rs = header->SectionNumber();
    if (rs < fMinSectionNumber)
      fMinSectionNumber = rs;
    if (ev < fMinEventNumber  )
      fMinEventNumber   = ev;
  }

  if (rn > fMaxRunNumber) {
    fMaxRunNumber   = rn;
    fMaxEventNumber = ev;
  }
  else if (rn == fMaxRunNumber) {
    if (ev > fMaxEventNumber)
      fMaxEventNumber = ev;
  }
//-----------------------------------------------------------------------------
//  run sections
//-----------------------------------------------------------------------------
  int ns = fCurrentRunRecord->fListOfRunSections.size();

  int found = 0;
  for (int i=0; i<ns; i++) {
    int rss = fCurrentRunRecord->fListOfRunSections[i];
    if (rs == rss) {
      found = 1;
      break;
    }
  }

  if (found == 0) {
    fCurrentRunRecord->fListOfRunSections.push_back(rs);
  }
  
  fNEvents++;

  return 0;
}


//_____________________________________________________________________________
int TDFCModule::EndJob() {

  char  filename[200];
  char  id[200];
  strncpy(id,fDatasetID.Data(),50);
  sprintf(filename,"%c%c%06x.%04x%1c%1c%1c%1c",
	  id[0],id[5],fMinRunNumber,fMinSectionNumber,id[1],id[2],id[3],id[4]);
  fNewFileName = TString(filename);

  //if (PrintLevel() > 0) 
  //  printf(" TDFCModule::EndJob - Number of events processed: %d\n",fNEvents);

  if (fNEvents <= 0)
    fFile = 0;

  if (PrintLevel() > 0)
    PrintCatalog(fPrintOpt);
 
  if (PrintLevel() > 10) {
    TString move, dfc;

    GetMvCommand(move);
    //    GetDfcCommand(dfc);
    //    printf("%s\n",move.Data());
  }

  if (fExecCommands)
    ExecCommands();

  return fReturnCode;
}

//_____________________________________________________________________________
int TDFCModule::GetNewFileName(TString& Name) {
  Name = fNewFileName;
  return 0;
}

//_____________________________________________________________________________
int TDFCModule::GetMvCommand(TString& Command) {
  // output dir may have form:
  // local : /a/b
  // remote: TUrl format, i.e. ftp://ewk@fcdfdata033.fnal.gov/gjet08/emoe
  // fNewFileName is always local to the directory

  TString dir;

  if (fOutputDir == NULL) {

					// local directory
    Command  = Form("mv %s %s",
		    fFileName.Data(),
		    fNewFileName.Data());
  }
  else if (strcmp(fOutputDir->GetHost(),"") == 0) {
					// local directory
    Command  = Form("mv %s  %s/%s ",
		    fFile->GetName(),
		    fOutputDir->GetFile(),
		    fNewFileName.Data());
  }
  else {
//-----------------------------------------------------------------------------
// final destination: remote server
// assume ftp works. also assume that the file is either 
// on a remote host or on a local one, but not on the 3rd one
//-----------------------------------------------------------------------------
    char  local_host[100];
    FILE* pipe = gSystem->OpenPipe("hostname -f","r");
    fscanf(pipe,"%s",local_host);
    gSystem->ClosePipe(pipe);
    
    TUrl file(fFile->GetName());

    if ((strcmp(file.GetHost(),"") == 0) || 
	 strcmp(local_host,file.GetHost()) == 0) {
//-----------------------------------------------------------------------------
//  file is local, rename, then scp 
//-----------------------------------------------------------------------------
      TString remote_copy;
      if ( gSystem->Getenv("USE_SCP") == 0)
	remote_copy = "fcp -c /usr/krb5/bin/rcp -N -r";
      else
	remote_copy = "scp";

      Command = Form("mv %s %s ; %s %s %s@%s:%s ; if [ $? == 0 ] ; then rm %s ; fi",
		      fFile->GetName(), 
		      fNewFileName.Data(),
                      remote_copy.Data(),
		      fNewFileName.Data(),
		      fOutputDir->GetUser(),
		      fOutputDir->GetHost(),
		      fOutputDir->GetFile(),
		      fNewFileName.Data());
    }
    else {
//-----------------------------------------------------------------------------
//  renamed file is located on a remote host, use ftp to rename it
//-----------------------------------------------------------------------------
      Command  = Form("ftp %s <<\r\n",fOutputDir->GetHost());
      Command += Form("%s\n",fOutputDir->GetUser());
      Command += Form("cd %s\n",fOutputDir->GetFile());
      Command += Form("rename %s %s\n",
		      gSystem->BaseName(fFile->GetName()),fNewFileName.Data());
    }

  }

  return 0;
}

//_____________________________________________________________________________
void TDFCModule::PrintFilename()
{
  // DatasetID is 6 characters long
  printf("new filename: %s \n",fNewFileName.Data());
  printf("run_min, event_min : %8i %8i \n",fMinRunNumber,fMinEventNumber);
  printf("run_max, event_max : %8i %8i \n",fMaxRunNumber,fMaxEventNumber);
  float fsize = fFile->GetSize()/1000000.;
  printf("file size          : %10.5f  \n",fsize);

  TString move, dfc;
  GetMvCommand(move);
  GetDfcCommand(dfc);
  printf("%s\n",move.Data());
}

//_____________________________________________________________________________
int TDFCModule::GetDfcCommand(TString& Cmd) {

  // DatasetID is 6 characters long
  float fsize = 0;
  if (fNewFile)
    fsize = fNewFile->GetSize()/1000000.;
  else 
    fsize = fFile   ->GetSize()/1000000.;

  Cmd  = Form("DFCFileTool -create -db %s -book %s ",
	      fDbID.Data(),fBook.Data());
  Cmd += Form("-file %s ",fNewFileName.Data()); 
  Cmd += Form("-dataset_id %s ",fDatasetID.Data());
  Cmd += Form("-id orphan -file_size %12.6f ",fsize);
  Cmd += Form("-event_count %i ",fNEvents);
  Cmd +=      "-lum_sum_online 0 -lum_sum_offline 0 ";
  Cmd += Form("-low_event \"[%i,%i]\" -high_event \"[%i,%i]\" ",
	      fMinRunNumber,fMinEventNumber,
	      fMaxRunNumber,fMaxEventNumber);
  Cmd += Form(" -runsecs \"\" -livetimes \"\" ");

  return 0;
}

//_____________________________________________________________________________
void TDFCModule::ExecCommands() {

  int i1/*, i2*/;
  
  TString move, dfc;

  if (fExecCommands) {
    i1  = (fExecCommands   ) % 10 ;
    //    i2  = (fExecCommands/10) % 10 ;
    
    if (i1 == 1) {
      GetMvCommand(move);
      printf(" **** executing move command: %s\n",move.Data());
      fReturnCode = gSystem->Exec(move.Data());
      if (fReturnCode == 0)
	fNewFile = new TFile(TString(TString(fOutputDir->GetFile())+
				     TString("/")+
				     fNewFileName.Data()));
    }

//     if ((fReturnCode == 0) && (i2 == 1)) {
//       GetDfcCommand(dfc);
//       printf(" **** executing DFC  command: %s\n",dfc.Data());
//       fReturnCode = gSystem->Exec(dfc.Data());
//     }
  }
}

//_____________________________________________________________________________
void TDFCModule::SetDataSet(const char* id, const char* book, const char* db) {
  fDatasetID = id;
  fBook      = book;
  fDbID      = db;
}

//_____________________________________________________________________________
void TDFCModule::PrintCatalog(const char* Opt) {
  // filename size date nevents low_run low_event high_run high_event
  char  cmd[1000], fn[1000];

  if ((strcmp(Opt,"") == 0) || (strstr(Opt,"banner") > 0)) {
    printf("----------------------------------------------------------------");
    printf("------------------------------------------------------------\n");
    printf("  MetaData    filename                         size       date  ");
    printf("  time   N(evts)   LowRun   LowEvent    HighRun  HighEvent  \n");
    printf("----------------------------------------------------------------");
    printf("------------------------------------------------------------\n");
  }

  if ((strcmp(Opt,"") == 0) || (strstr(Opt,"data") > 0)) {

    sprintf(cmd,"echo %s | awk -F / '{print $NF}'",fFileName.Data());
    FILE* f = gSystem->OpenPipe(cmd,"r");
    fscanf(f,"%s",fn);
    fclose(f);

    //if (fFile)
      printf("%s.%04d %-20s %10.3f %20s %6i %7i %9i %7i %9i\n",
	     fDatasetID.Data(),
	     fNFileset,
	     fn,
	     fFileSize,
	     fFileDate.Data(),
	     //fFile->GetSize()/1000000.,
	     //fFile->GetModificationDate().AsSQLString(),
	     fNEvents,
	     fMinRunNumber,
	     fMinEventNumber,
	     fMaxRunNumber,
	     fMaxEventNumber);
//    else
//      printf("%s.%04d %-20s %10.3f %20s %6i %7i %9i %7i %9i\n",
//	     fDatasetID.Data(),
//	     fNFileset,
//	     fn,
//	     0.,                     // making up bogus numbers
//	     "1970-01-01 00:00:00",  // making up bogus numbers
//	     fNEvents,
//	     fMinRunNumber,
//	     fMinEventNumber,
//	     fMaxRunNumber,
//	     fMaxEventNumber);
//      
  }

  int nr     = fRunRecord.size();
  if (PrintLevel() % 10 == 2) {
//-----------------------------------------------------------------------------
//  print run-level catalog
//-----------------------------------------------------------------------------
    int offset = fDatasetID.Length()+21+11+21+5;
    for (int i=0; i<nr; i++) {
      RunRecord_t *r = fRunRecord[i];
      for (int i=0; i<offset; i++) printf(" ");
      printf(" %6i %7i %9i %7i %9i\n",
	     r->fNEvents,r->fNumber,r->fLoEvent,r->fNumber,r->fHiEvent);
    }
    
    if (PrintLevel() != 0 &&
	((strcmp(Opt,"") == 0) || (strstr(Opt,"foot") > 0))) {
      printf("--------------------------------------------------------------");
      printf("----------------------------------------------------------\n");
    }
  }

  else if (PrintLevel() % 10 == 3) {
//-----------------------------------------------------------------------------
//  print getFileInfo style catalog assuming that we run on a single file
//-----------------------------------------------------------------------------
    printf("<getfileinfo>\n");
    printf("name           : %-s\n",fFileName.Data());
    printf("size           : %lld \n",fFile->GetSize());
    printf("First run/event: %6i / %7i \n", fMinRunNumber,fMinEventNumber);
    printf("Last  run/event: %6i / %7i \n", fMaxRunNumber,fMaxEventNumber);
    printf("events         : %6i \n",fNEvents);
    printf("runsections    : ");

    for (int i=0; i<nr; i++) {
      RunRecord_t *r = fRunRecord[i];
      //      printf("run number: %i\n",r->fNumber);
      int ns = r->fListOfRunSections.size();
//-----------------------------------------------------------------------------
// order the runsections
//-----------------------------------------------------------------------------
      int* a    = new int[ns];
      int* rr   = new int[2*ns];

      for (int i=0; i<ns; i++) {
	a[i] = r->fListOfRunSections[i];
	for (int j=i+1; j<ns; j++) {
	  int aj = r->fListOfRunSections[j];
	  if (aj < a[i]) {
	    r->fListOfRunSections[j] = a[i];
	    a[i]   = aj;
	  }
	}
	//	printf("sec number: %i\n",a[i]);
      }
//-----------------------------------------------------------------------------
// run sections ordered, make ranges
//-----------------------------------------------------------------------------
      int nr    = 0;
      
      for (int i=0; i<ns; i++) {
	if ((nr > 0) && (a[i] - rr[2*nr-1] == 1)) {
				// extend the existing range
	  rr[2*nr-1] = a[i];
	}
	else {
				// new range
	  nr        += 1;
	  rr[2*nr-2] = a[i];
	  rr[2*nr-1] = a[i];
	}
      }
//-----------------------------------------------------------------------------
//  print run section ranges
//-----------------------------------------------------------------------------
      for (int i=0; i<nr; i++) {
	printf("%6i/%i:%i ",r->fNumber,rr[2*i],rr[2*i+1]);
      }
      delete [] a;
      delete [] rr;
    }
    printf("\n");
    printf("</getfileinfo>\n");
    
    if (PrintLevel() != 0 &&
	((strcmp(Opt,"") == 0) || (strstr(Opt,"foot") > 0))) {
      printf("--------------------------------------------------------------");
      printf("----------------------------------------------------------\n");
    }
  }
}
