//
#ifndef CalPatRec_HlPrint_hh
#define CalPatRec_HlPrint_hh

#include "TObject.h"
#include "TObjArray.h"
#include "TString.h"

#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"

#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"

namespace mu2e {

  class StrawHit;
  class CaloCluster;
  class CaloProtoCluster;
  class CrvDigi;
  class CrvRecoPulse;
  class CrvCoincidence;
  class CrvCoincidenceCluster;
  class GenParticle;
  class TimeCluster;
  class KalSeed;
  class ComboHit;
  class HelixSeed;
  class SimParticle;
  class StepPointMC;
  class StrawDigiMC;
  class StrawGasStep;
  class TrackClusterMatch;
  class TrkCaloHit;
  class TrkStrawHit;
  class TrkPrintUtils;

class HlPrint {
public:

  const art::Event*   _event;
  double              _tmp[100];  // for testing
  TrkPrintUtils*      _printUtils;

private:

  HlPrint(const fhicl::ParameterSet* Pset = NULL);
  ~HlPrint();

  class  Cleaner {
  public:
    Cleaner();
    ~Cleaner();
  };

  friend class Cleaner;

  static HlPrint*  _Instance;
public:
//-----------------------------------------------------------------------------
// HlPrint
//-----------------------------------------------------------------------------
  static HlPrint*  Instance(const fhicl::ParameterSet* PSet = NULL);
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
  const art::Event*              Event      () { return _event      ; }
//-----------------------------------------------------------------------------
// other methods
//-----------------------------------------------------------------------------
  void   AddObject      (const char* Name, void* Object);
  void*  FindNamedObject(const char* Name);

  void   SetEvent(const art::Event* Evt) { _event = Evt; }

  void   printEventHeader();
//-----------------------------------------------------------------------------
// calorimeter
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// CRV
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// tracking
//-----------------------------------------------------------------------------
  void printComboHit      (const ComboHit*     Hit       ,
                           const StrawGasStep* Step      ,
                           const char*         Opt   = "",
                           int                 INit  = -1,
                           int                 Flags = -1);

  void printComboHitCollection (const char* ChCollTag              ,  // ComboHit     collection
                                const char* ChfCollTag             ,  // StrawHitFlag collection
                                const char* SdmcMCCollTag = nullptr,  // "makeSD" or "compressDigiMCs"
                                double      TMin          = -1.e6  ,
                                double      TMax          =  1.e6  );
//-----------------------------------------------------------------------------
// MC truth: gen and sim particles
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// pass the detector name to know what to print for different detectors
// tested for Detector = 'tracker', 'calorimeter'
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// time clusters
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// calorimeter cluster added to the track fit
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// extrapolation and track-to-calorimeter matching
//-----------------------------------------------------------------------------
};
}

#endif
