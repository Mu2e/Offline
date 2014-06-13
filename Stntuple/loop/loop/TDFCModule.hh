#ifndef TDFCModule_hh
#define TDFCModule_hh

#include <vector>
#include "TUrl.h"
#include "Stntuple/loop/TStnModule.hh"

class TFile;

class TDFCModule: public TStnModule {
  // everything is public, it is communism
public:

  struct RunRecord_t {
    Int_t       fNumber;                // run number
    Int_t       fNEvents;               // number of events
    Int_t       fLoEvent;	        // low event number
    Int_t       fHiEvent;	        // high event number
    std::vector<int> fListOfRunSections;     // run section ranges
  };

  Int_t         fMinRunNumber;
  Int_t         fMinSectionNumber;
  Int_t         fMinEventNumber;
  Int_t         fMaxRunNumber;
  Int_t         fMaxEventNumber;
  Int_t         fNEvents;
  TString       fDatasetID;
  TString       fBook;
  TString       fDbID;
  TString       fPrintOpt;
  Int_t         fNFileset;
  Int_t         fExecCommands;
  TUrl*         fOutputDir;     //! if move to other directory(ftp, if remote)
  RunRecord_t*  fCurrentRunRecord;
  std::vector<RunRecord_t*>
                fRunRecord;
  Int_t         fReturnCode;

private:
  TFile*        fFile;
  TString       fFileName;
  TString       fFileDate;
  Float_t       fFileSize;
  TFile*        fNewFile;
  TString       fNewFileName;

//-----------------------------------------------------------------------------
//  functions
//-----------------------------------------------------------------------------
public:
  TDFCModule(const char* name = "DFC", const char* title = "DFC Module");
  ~TDFCModule();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  int       GetNewFileName (TString& Name   );
  int       GetMvCommand   (TString& Command);
  int       GetDfcCommand  (TString& Command);

  int       ReturnCode     () { return fReturnCode; }
  TString   GetFileDate    () { return fFileDate;}
  TString   GetNewFileName () { return fNewFileName;}
  TString   GetOldFileName () { return fFileName;}
  Float_t   GetFileSize    () { return fFileSize;} // MB
//-----------------------------------------------------------------------------
// setters ... SetOutputDir leaks memory. I know.
//-----------------------------------------------------------------------------
  void      SetExecCommands(Int_t Flag) { fExecCommands = Flag; }
  void      SetDataSet     (const char* DsID, const char* Book, const char* Db);
  void      SetOutputDir   (const char* Dir) { fOutputDir = new TUrl(Dir) ; }
  void      SetPrintOpt    (const char* Opt) { fPrintOpt = Opt; }
  void      SetNFileset    (Int_t n) { fNFileset = n; }
//-----------------------------------------------------------------------------
// other methods
//-----------------------------------------------------------------------------
  void      PrintFilename  ();
  void      PrintCatalog   (const char* Opt="");
  void      ExecCommands   ();
//-----------------------------------------------------------------------------
// overloaded methods of TStnModule
//-----------------------------------------------------------------------------
  int       BeginJob       ();
  int       BeginRun       ();
  int       Event          (int ientry);
  int       EndJob         ();

  ClassDef(TDFCModule,0)
};
#endif
