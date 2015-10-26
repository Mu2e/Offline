///////////////////////////////////////////////////////////////////////////////
//
//

#include "Stntuple/mod/TModule.hh"
#include "Stntuple/mod/TAnaRint.hh"
#include "Stntuple/mod/TAnaDump.hh"


// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

// C++ includes.
#include <iostream>

#include "TString.h"

//-----------------------------------------------------------------------------
TModule::TModule(fhicl::ParameterSet const& PSet, const char* Name):
  TNamed(Name,Name)
{

  int    n, index;
  static char* dummy[100];

  fFile   = 0;
  fFolder = new TFolder(GetName(),GetName());

  memset(fDebugBit,0,kNDebugBits*sizeof(int));

  fFclDebugBits    = PSet.get<fhicl::ParameterSet>("debugBits"                );
  fInteractiveMode = PSet.get<int>                ("interactiveMode",        0);

  fAnaRint         = TAnaRint::Instance(0,dummy);
  fDump            = TAnaDump::Instance();

  const char* key;
                                        // a flag is an integer!
  n = fFclDebugBits.get_names().size();
  for (int i=0; i<n; i++) {
    key                = fFclDebugBits.get_names().at(i).data();
    sscanf(key,"bit%i" ,&index );

    fDebugBit[index]  = fFclDebugBits.get<int>(key);

    printf("... TModule: bit=%3i is set to %i \n",index,fDebugBit[index]);
  }

};

//-----------------------------------------------------------------------------
TModule::~TModule() {
  if (fFile) {
    delete fFile;
    fFile  = NULL;
  }
};


//-----------------------------------------------------------------------------
void TModule::beginJob() {
};


//-----------------------------------------------------------------------------
bool TModule::beginRun(art::Run &  Rn) {
  return true;
};

//-----------------------------------------------------------------------------
bool TModule::filter  (art::Event& Evt) {
  if (fInteractiveMode != 0) {
    fDump->SetEvent(Evt);
    fAnaRint->SetInteractiveMode(fInteractiveMode);
    fAnaRint->Rint()->Run(true);
                                        // provision for switching the interactive mode OFF

    fAnaRint->GetInteractiveMode(fInteractiveMode);
  }
  return true;
};

//-----------------------------------------------------------------------------
void TModule::endJob  () {
};

//_____________________________________________________________________________
  void     TModule::AddHistogram(TObject* hist, const char* FolderName) {
    TFolder* fol = (TFolder*) fFolder->FindObject(FolderName);
    fol->Add(hist);
  }

//_____________________________________________________________________________
  void TModule::HBook1F(TH1F*& Hist, const char* Name, const char* Title,
                           Int_t Nx, Double_t XMin, Double_t XMax,
                           const char* FolderName)
  {
    // book 2D histogram, add it to the module's list of histograms and
    // return pointer to it to the user

    Hist = new TH1F(Name,Title,Nx,XMin,XMax);
    AddHistogram(Hist,FolderName);
  }

//_____________________________________________________________________________
  void TModule::HBook2F(TH2F*& Hist, const char* Name, const char* Title,
                           Int_t Nx, Double_t XMin, Double_t XMax,
                           Int_t Ny, Double_t YMin, Double_t YMax,
                           const char* FolderName)
  {
    // book 2D histogram, add it to the module's list of histograms and
    // return pointer to it to the user

    Hist = new TH2F(Name,Title,Nx,XMin,XMax,Ny,YMin,YMax);
    AddHistogram(Hist,FolderName);
  }

//_____________________________________________________________________________
  void TModule::HProf(TProfile*& Hist, const char* Name, const char* Title,
                         Int_t Nx, Double_t XMin, Double_t XMax,
                         Double_t YMin, Double_t YMax,
                         const char* FolderName)
  {
    // book 2D histogram, add it to the module's list of histograms and
    // return pointer to it to the user

    Hist = new TProfile(Name,Title,Nx,XMin,XMax,YMin,YMax);
    AddHistogram(Hist,FolderName);
  }


  //_____________________________________________________________________________
  int  TModule::SaveFolder(TFolder* Folder, TDirectory* Dir) {
  // save Folder into a subdirectory
  // do not write TStnModule's - for each TStnModule save contents of its
  // fFolder

    //    TFolder*     fol;
    TDirectory*  dir;
    TObject*     o;
//-----------------------------------------------------------------------------
// create new subdirectory in Dir to save Folder
//-----------------------------------------------------------------------------
    Dir->cd();
    dir = Dir->mkdir(Folder->GetName(),Folder->GetName());
    dir->cd();

    TIter  it(Folder->GetListOfFolders());
    while ((o = it.Next())) {
      if (strcmp(o->ClassName(),"TFolder") == 0) {
        SaveFolder((TFolder*) o, dir);
      }
      else if (! o->InheritsFrom("TStnModule")) {
        o->Write();
        //      gDirectory->GetListOfKeys()->Print();
      }
    }

    Dir->cd();
    return 0;
  }
