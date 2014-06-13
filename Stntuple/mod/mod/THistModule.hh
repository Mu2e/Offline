//--------------------------------------------------------------------------
// File and Version Information: RootHistModule.hh,v 1.0
//
// Description:
//	Class THistModule: base class for ROOT histogramming modules
//
// Environment: CDF Run II
//
// Author List: P.Murat
//
// Copyright Information: 
//   Copyright (C) 1999		CDF/Fermilab
//------------------------------------------------------------------------

#ifndef THistModule_HH
#define THistModule_HH

#include "TString.h"
#include "TFile.h"
#include "TObjArray.h"

#include "Stntuple/obj/AbsEvent.hh"
#include "Stntuple/mod/TModule.hh"

class TTree;

class THistModule : public TModule {
//------------------------------------------------------------------------------
//  static data members
//------------------------------------------------------------------------------
protected:
					// there are some initializations 
					// which need to be done just once
  static TString    fgFileName;
  static TFile*     fgFile;
  static TObjArray* fgModuleList;
  static TTree*     fgTree;
  static int        fgMakeSubdirs;
  static int        fgMaxFileSize;
  static int        fgFileNumber;
  static int        fgOpenNextFile;
  static int        fgSplitLevel;
  static int        fgBufferSize;
  static int        fgCompressionLevel;
//------------------------------------------------------------------------------
//  data members of the module
//------------------------------------------------------------------------------
					// name of the directory in a ROOT file
					// associated with the module
  TString       fDirName;
					// list of histograms/ntuples owned by
					// the module
  TObjArray*    fHistogramList;
					// cache for cd() command
  TDirectory*   fOldDir;
//------------------------------------------------------------------------------
//  methods of the class
//------------------------------------------------------------------------------
public:
					// ****** constructors and destructor

  THistModule(fhicl::ParameterSet const& PSet, const char* Name);

  ~THistModule( );
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  int        SplitLevel      () { return fgSplitLevel;   }
  int        CompressionLevel() { return fgCompressionLevel; }
  int        BufferSize      () { return fgBufferSize;       }
  TFile*     File            () { return fgFile;         }
  TObjArray* HistogramList   () { return fHistogramList; }

//-----------------------------------------------------------------------------
// other methods
// define name of the output ROOT file, this can be done once per job
//-----------------------------------------------------------------------------
  int        SetFileName    (const char* Filename);
  int        OpenNewFile    (const char* Filename);

					// ****** overloaded methods of the 
					// base class

					// use beginJob to book histograms
					// use event to fill histograms

					// ****** methods called by the action
					// controller

  virtual int beforeBeginJob();
  virtual int afterBeginJob ();
  virtual int beforeBeginRun(art::Run& aRun );
  virtual int afterBeginRun (art::Run& aRun );
  virtual int beforeEvent   (AbsEvent& event);
  virtual int afterEvent    (AbsEvent& event);
  virtual int beforeEndRun  (art::Run& aRun );
  virtual int afterEndRun   (art::Run& aRun );
  virtual int beforeEndJob  ();
  virtual int afterEndJob   ();

  //  ClassDef(THistModule,0)
};

#endif
