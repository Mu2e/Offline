///////////////////////////////////////////////////////////////////////////////
//
//

// Mu2e includes.
// #include "RecoDataProducts/inc/StrawHitCollection.hh"
// #include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
// #include "RecoDataProducts/inc/CaloHitCollection.hh"
// #include "RecoDataProducts/inc/CaloClusterCollection.hh"
// #include "MCDataProducts/inc/GenParticleCollection.hh"

// Framework includes.
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"

// C++ includes.
#include <iostream>

#include "TString.h"
#include "TFolder.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

// namespace murat {


class TAnaRint;
class TAnaDump;

class TModule : public art::EDFilter, public TNamed {

  enum { kNDebugBits = 100 };

public:
				// there are some initializations which need 
				// to be done just once

  static int          fgInitialized;

  TFile*              fFile;

  int                 fDebugBit[kNDebugBits];		// flags for different debug options

  fhicl::ParameterSet fFclDebugBits;
  int                 fInteractiveMode;

  int                 fEventNumber;
  int                 fRunNumber;
  int                 fPassed;                          // for filtering

  TAnaRint*           fAnaRint;
  TAnaDump*           fDump;

  art::Run*           fRun;
				// provides for a possibility for any ROOT 
				// module to call (whenever necessary) a
				// loaded in interpreted function which 
				// name has to be defined at run time
				// it is a user responsibility to load a
				// macro with this function. Its signature:
				// 
				//       int function_name(int flag)
				// 
				// fFunctionName.Data() would return the 
				// `function_name'
  TString    fFunctionName;
				// each TModule adds a folder to gROOT
  TFolder*            fFolder;
//-----------------------------------------------------------------------------
// methods - do not overload ::filter...
//-----------------------------------------------------------------------------
  explicit TModule(fhicl::ParameterSet const& pset, const char* Name);
  virtual ~TModule();

  virtual void beginJob();
  virtual bool beginRun(art::Run &  Rn);
  virtual bool filter  (art::Event& Evt);
  virtual void endJob  ();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  TFolder* GetFolder() { return fFolder; }

  int      DebugBit(int I) { return fDebugBit[I]; }

  const char* FunctionName() { return fFunctionName.Data(); }
//-----------------------------------------------------------------------------
// the following helper methods allow to save 1 line per request, which in 
// case of 100's histograms booked is a non-negligible number
//-----------------------------------------------------------------------------
  void  AddHistogram(TObject* hist, const char* FolderName = "Hist");
  
  void  HBook1F(TH1F*& Hist, const char* Name, const char* Title,
		Int_t Nx, Double_t XMin, Double_t XMax,
		const char* FolderName = "Hist");
  
  void  HBook2F(TH2F*& Hist, const char* Name, const char* Title,
		Int_t Nx, Double_t XMin, Double_t XMax,
		Int_t Ny, Double_t YMin, Double_t YMax,
		const char* FolderName = "Hist");

  void  HProf (TProfile*& Hist, const char* Name, const char* Title,
	       Int_t Nx, Double_t XMin, Double_t XMax,
	       Double_t YMin, Double_t YMax,
	       const char* FolderName = "Hist");
  
  int  SaveFolder(TFolder* Folder, TDirectory* Dir);
  
};

// }  // end namespace murat
