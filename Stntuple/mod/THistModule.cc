//--------------------------------------------------------------------------
// Description:
//	Class THistModule: base class for ROOT histogramming modules
//
// Environment: Mu2e
//
// Author List: P.Murat
//
// Copyright Information: 
//   Copyright (C) 1999		CDF/Fermilab
//------------------------------------------------------------------------

#include <assert.h>
#include <string>

#include "TInterpreter.h"
#include "TH1.h"
#include "TFolder.h"

#include "Stntuple/obj/AbsEvent.hh"
#include "Stntuple/mod/THistModule.hh"

// ClassImp(THistModule)
					// initialize static class members

TObjArray* THistModule::fgModuleList   = 0;
TFile*     THistModule::fgFile         = 0;
TTree*     THistModule::fgTree         = 0;
int        THistModule::fgMakeSubdirs  = 1;

int        THistModule::fgBufferSize   = 64000;  // in bytes
int        THistModule::fgMaxFileSize  =   350;  // in MBytes
TString    THistModule::fgFileName     = "";
int        THistModule::fgSplitLevel   = 99;
int        THistModule::fgCompressionLevel = 1;

int        THistModule::fgFileNumber   = 0;
int        THistModule::fgOpenNextFile = 0;


//______________________________________________________________________________
THistModule::THistModule(fhicl::ParameterSet const& PSet, const char* Name): 
  TModule(PSet,Name)
{
  fOldDir        = NULL;
  fHistogramList = new TObjArray(10);

  fgBufferSize       = PSet.get<int>        ("bufferSize"      ,fgBufferSize);
  fgMaxFileSize      = PSet.get<int>        ("maxFileSize"     ,fgMaxFileSize);
  fgFileName         = PSet.get<std::string>("histFileName"    ,fgFileName.Data()).data();
  fgSplitLevel       = PSet.get<int>        ("splitLevel"      ,fgSplitLevel);
  fgCompressionLevel = PSet.get<int>        ("compressionLevel",fgCompressionLevel);

  fHistogramList->SetName("HistogramList");
  fFolder->Add(fHistogramList);

  fFolder->AddFolder("Hist","Hist");

  if (fgModuleList == 0) {
    fgModuleList = new TObjArray(10);
    //    framework()->actions()->append(new THistModuleAction());
//-----------------------------------------------------------------------------
// do not append the histograms to the directory list of objects to allow
// using short histogram names
//-----------------------------------------------------------------------------
    TH1::AddDirectory(kFALSE);
  }
  fgModuleList->Add(this);
}


//______________________________________________________________________________
THistModule::~THistModule() {

//-----------------------------------------------------------------------------
//  close the file in the destructor such that all THistModules would
//  have their beforeEndJob entry points executed
//-----------------------------------------------------------------------------

  if (fgFile) {
    fgFile->Write();
    delete fgFile;
    fgFile = 0;
  }
  if (fHistogramList) {
    delete fHistogramList;
    fHistogramList = 0;
  }

  if (fgModuleList) {
    delete fgModuleList;
    fgModuleList = 0;
  }
}


//______________________________________________________________________________
int THistModule::OpenNewFile(const char* Filename) {
  // open next file
  int rc = 0;
  fgFile = new TFile(Filename,"RECREATE");
  if (! fgFile) {
    Error("beginJob","an attempxt to open a new ROOT file %s failed",
	  fgFileName.Data());
    rc = -1;
  }
  return rc;
}


//______________________________________________________________________________
int THistModule::beforeBeginJob() {

					// return code
  int rc  =  0;
				// give more time to define TModule::fName
  fDirName = GetName();

  if ((fgFileName != "") && (fgFile == 0)) {

					// create a new ROOT File 
    rc = OpenNewFile(fgFileName.Data());
  }

  if (fgFile) {
					// file opened, don't forget to make
					// new directory before booking 
    if (fgMakeSubdirs) {
      fOldDir = gDirectory;
      fgFile->mkdir(fDirName);
      fgFile->cd(fDirName);
    }
  }
  return rc;
}


//______________________________________________________________________________
int THistModule::afterBeginJob() {

  char par[200];

  if (strcmp(FunctionName(),"") != 0) {
					// call interpreted function if it 
					// is defined
    sprintf(par,"0,0x%lx",(long int) this);
    gInterpreter->Execute(FunctionName(),par);
  }

  if (fgMakeSubdirs && fgFile) {
    fOldDir->cd();
  }

  return 0;
}


//______________________________________________________________________________
int THistModule::beforeBeginRun(art::Run& aRun) {
  return 0;
}
    

//______________________________________________________________________________
int THistModule::afterBeginRun(art::Run& aRun) {
  return 0;
}
    

//______________________________________________________________________________
int THistModule::beforeEvent(AbsEvent& event) {
  if (fgFile) {
					// don't forget to change the directory
					// before starting filling histograms
    fOldDir = gDirectory;

    if (fgOpenNextFile) {
//-----------------------------------------------------------------------------
// current file will be closed by FillStntupleModule, this event will be 
// written into the next one, so write out the folder, 
// then reset the histograms
//-----------------------------------------------------------------------------
      fgFile->cd();
      fFolder->Write();
      fOldDir->cd();

      TObjArray* list = (TObjArray*) fFolder->FindObject("Hist");
      if (list) {
	int nhist = list->GetEntriesFast();
	for (int i=0; i<nhist; i++) {
	  TH1* h = (TH1*) list->UncheckedAt(i);
	  h->Reset();
	}
      }
    }
//-----------------------------------------------------------------------------
//  I suspect fgMakeSubdirs is always 1 now
//-----------------------------------------------------------------------------
    if (fgMakeSubdirs) {
      fgFile->cd(fDirName);
    }
  }
  return 0;
}
    

//______________________________________________________________________________
int THistModule::afterEvent(AbsEvent& event) {
  // need to cd() back 

  char par[200] ;
  if (strcmp(FunctionName(),"") != 0) {
					// call interpreted function if it 
					// is defined

    sprintf(par,"1,0x%lx",(long int) this);
    gInterpreter->Execute(FunctionName(),par);
  }

  if (fgFile) {
    if (fgMakeSubdirs) {
      fOldDir->cd();
    }
  }
  return 0;
}
    

//______________________________________________________________________________
int THistModule::beforeEndRun(art::Run& aRun) {
  return 0;
}
    

//______________________________________________________________________________
int THistModule::afterEndRun(art::Run& aRun) {
  return 0;
}
    

//______________________________________________________________________________
int THistModule::beforeEndJob() {
  if (fgFile) {
					// don't forget to change the directory
					// before starting filling histograms
    fOldDir = gDirectory;

    fgFile->cd();
    fFolder->Write();
    fOldDir->cd();
  }

  return 0;
}
    

//-----------------------------------------------------------------------------
int THistModule::afterEndJob() {
					// return code
  int   rc = 0;
  char  par[200];

  if (strcmp(FunctionName(),"") != 0) {
					// call interpreted function if it is defined
    sprintf(par,"2,0x%lx",(long int) this);
    gInterpreter->Execute(FunctionName(),par);
  }

  return rc;
}
    

//_____________________________________________________________________________
int THistModule::SetFileName(const char* nm) {
  if (fgFile) {
    printf(">>> ERROR: ROOT histogram file is already open. Renaming is disabled.\n");
  }
  else {
    if (nm && (fgFileName  != "")) {
      printf(">>> WARNING : %s is setting a new name for ROOT histogram file : %s\n",
	     GetName(), nm);
    }
    fgFileName = nm;
  }
  return 0;
}

