//-----------------------------------------------------------------------------
//  Dec 28 2004 P.Murat: base class for STNTUPLE input module
//-----------------------------------------------------------------------------
#include "TChain.h"
#include "TFile.h"
#include "TSystem.h"
#include "TRegexp.h"
#include "TObjString.h"

#include "Stntuple/base/TStnDataset.hh"
#include "Stntuple/base/TStnFileset.hh"
#include "Stntuple/base/TCdf2Files.hh"

ClassImp(TStnDataset)

//_____________________________________________________________________________
TStnDataset::TStnDataset(const char* Name): TNamed(Name,Name) {
  // non-cataloged dataset

  fNEvents        = 0;
  fMinRunNumber   = -1;
  fMaxRunNumber   = 10000000;
  fBook           = "";
  fChain          = 0;
  fListOfFiles    = new TObjArray();
  fListOfFilesets = new TObjArray();
  fListOfBadFiles = new TObjArray();
  fCataloged      = 0;
  fMcFlag         = -1;
  fFilesetFormat  = "";
  fDoneBadFiles   = 0;
}

//_____________________________________________________________________________
TStnDataset::TStnDataset(const char* Book,
			 const char* Name  , 
			 Int_t       MinRunNumber, 
			 Int_t       MaxRunNumber,
			 const char* Type        ):
  TNamed(Name,Name)
{
  // cataloged dataset:Name,Book as it comes is like "stntuple/dev_240",gmbs08"

  fListOfFiles    = new TObjArray();
  fListOfFilesets = new TObjArray();
  fListOfBadFiles = new TObjArray();
  Init(Book,Name,MinRunNumber,MaxRunNumber,Type);
}

//_____________________________________________________________________________
TStnDataset::TStnDataset(const char* Book   ,
			 const char* Name   , 
			 const char* Fileset, 
			 const char* File   ,
			 const char* Type   ): 
  TNamed(Name,Name)
{
  // cataloged dataset:Name,Book as it comes is like "stntuple/dev_240",gmbs08"

  TString type = Type;

  fNEvents  = 0;
  fMinRunNumber = 1;
  fMaxRunNumber = 10000000;
  fBook         = Book;
  fChain        = 0;

  if (type == "STNTUPLE") {
    fChain    = new TChain("STNTUPLE");
  }

  fListOfFiles    = new TObjArray();
  fListOfFilesets = new TObjArray();
  fListOfBadFiles = new TObjArray();
  fDoneBadFiles   = 0;

  if (strcmp(Fileset,"") != 0) {
//-----------------------------------------------------------------------------
// so far handle only one fileset, in  principle this can be a list of fileset 
// names separated by ":" server=0 marks fileset as uninitialized
//-----------------------------------------------------------------------------
    TStnFileset* fs = new TStnFileset(Fileset,0,0,0,-1,-1);
    fListOfFilesets->Add((TObject*) fs);
  }

  fCataloged = 1;
}


//_____________________________________________________________________________
TStnDataset::~TStnDataset() {
  // it is more efficient to delete objects in the order opposite to 
  // their creation

  fListOfFiles->Delete();
  delete fListOfFiles;

  fListOfFilesets->Delete();
  delete fListOfFilesets;

  fListOfBadFiles->Delete();
  delete fListOfBadFiles;

  if (fChain   ) delete fChain   ;
}


//_____________________________________________________________________________
Int_t TStnDataset::Init(const char* Book,
			const char* Name  , 
			Int_t       MinRunNumber, 
			Int_t       MaxRunNumber,
			const char* Type        )
{
  // cataloged dataset:Name,Book as it comes is like "stntuple/dev_240",gmbs08"

  fName         = Name;
  fTitle        = Name;
  fNEvents      = 0;
  fMinRunNumber = MinRunNumber;
  fMaxRunNumber = MaxRunNumber;
  fBook         = Book;
  fChain        = 0;

  TString type = Type;
  if (type == "STNTUPLE") {
    fChain    = new TChain("STNTUPLE");
  }

  fDoneBadFiles   = 0;

  if ((strcmp(Book,"file") == 0) || (strcmp(Book,"dir") == 0)) fCataloged = 0;
  else                                                         fCataloged = 1;

  return 0;
}



//-----------------------------------------------------------------------------
Int_t TStnDataset::AddFileset(const char* Name) {

  TStnFileset* fs;

  TObject* found = fListOfFilesets->FindObject(Name);

  if (! found) {
    fs = new TStnFileset(Name,0,0,0,-1,-1);
    fListOfFilesets->Add((TObject*) fs);
  }

  return 0;
}

//-----------------------------------------------------------------------------
Int_t TStnDataset::AddFile(const char* Name) {
  // this method should be called only for non-cataloged datasets
  // make sure we're not adding the same file twice
  // as such, assume that the file is OK
  
  int nev;

  if (fCataloged == 1) {
    Error("AddFile(name)","can only be used for non-cataloged datasets");
    return -1;
  }

  TObject* found = fChain->GetListOfFiles()->FindObject(Name);

  if (! found) {
					// create new piece of metadata
    TFile* f    = TFile::Open(Name);
    TTree* tree = (TTree*) f->Get("STNTUPLE");
    nev         = int(tree->GetEntries());
    fNEvents   += nev;

    TCdf2Files* file = new TCdf2Files();

    file->fFILE_NAME    = Name;
    file->fFILESET_NAME = "none";
    file->fFILE_SIZE    = f->GetSize();
    file->fEVENT_COUNT  = nev;
    file->fLOW_EVENT    = 0;
    file->fHIGH_EVENT   = 100000000;
    file->fLOW_RUN      = 0;
    file->fHIGH_RUN     = 100000000;
    file->fSTATUS       = 0;

    fListOfFiles->Add((TObject*) file);
    
    fChain->AddFile(Name,nev);
    f->Close();
  }

  return 0;
}


//-----------------------------------------------------------------------------
Int_t TStnDataset::AddFile(const char* Name, 
			   const char* FilesetName,
			   Float_t     Size       ,
			   Int_t       NEvents    ,
			   Int_t       LoEvt      ,
			   Int_t       LoRun      ,
			   Int_t       HiEvt      ,
			   Int_t       HiRun      ,
			   Int_t       StatusCode)
{

  // this method should be called only for cataloged datasets

  if (fCataloged == 0) {
    Error("AddFile(name,fset..)","can only be used for cataloged datasets");
    return -1;
  }
					// create new piece of metadata
  TCdf2Files* file = new TCdf2Files();

  file->fFILE_NAME    = Name;
  file->fFILESET_NAME = FilesetName;
  file->fFILE_SIZE    = Size;
  file->fEVENT_COUNT  = NEvents;
  file->fLOW_EVENT    = LoEvt;
  file->fHIGH_EVENT   = HiEvt;
  file->fLOW_RUN      = LoRun;
  file->fHIGH_RUN     = HiRun;
  file->fSTATUS       = StatusCode;

  fListOfFiles->Add((TObject*) file);
//-----------------------------------------------------------------------------
// fNEvents - number of events in the chain!
//-----------------------------------------------------------------------------
  if (file->fSTATUS >= 0) {
    fChain->AddFile(Name,NEvents);
    fNEvents  += NEvents;
  }


  return 0;
}

//-----------------------------------------------------------------------------
Int_t TStnDataset::AddBadFile(const char* Name) {
  GetListOfBadFiles()->Add(new TObjString(Name));
  return 0;
}

//-----------------------------------------------------------------------------
void TStnDataset::Clear(Option_t* Opt) {
}

//_____________________________________________________________________________
void TStnDataset::Print(Option_t* Opt) const {

  TCdf2Files*  f;

  int nfiles = fListOfFiles->GetEntriesFast();
  
  for (int i=0; i<nfiles; i++) {
    f = (TCdf2Files*) fListOfFiles->At(i);
    if (i == 0) f->Print("banner");
    f->Print("data");
  }
}

