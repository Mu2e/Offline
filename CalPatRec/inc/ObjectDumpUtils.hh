//
#ifndef __CalPatRec_inc_ObjectDumpUtils_hh__
#define __CalPatRec_inc_ObjectDumpUtils_hh__

#include "TObject.h"
#include "TObjArray.h"
#include "TString.h"
#include "TGraph.h"
#include "TMarker.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"

namespace art {
  class Event;
}

class KalRep;

#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "RecoDataProducts/inc/CaloProtoClusterCollection.hh"

namespace mu2e {
  class StrawHit;
//  class StrawHitPosition;
  class CaloCluster;
  class CaloProtoCluster;
  class StepPointMC;
  class GenParticle;
  class SimParticle;
  //  class CalTimePeak;
  class TrkStrawHit;


  class ObjectDumpUtils {
  protected:

    static std::string                            _FlagBgrHitsModuleLabel;
    static const  StrawDigiMCCollection*          _ListOfMCStrawHits;
    static const  SimParticleTimeOffset*          _TimeOffsets;

  public:
//-----------------------------------------------------------------------------
// methods
//-----------------------------------------------------------------------------

    static void SetFlagBgrHitsModuleLabel(const char* Label) { _FlagBgrHitsModuleLabel = Label; }

    static void printEventHeader(const art::Event* Event, const char* Message = "");

    static void printKalRep(const KalRep* Krep, const char* Opt = "", const char* Prefix = "");

    static void printKalRepCollection(const art::Event*          Event        , 
				      const KalRepPtrCollection* Coll         ,
				      int                        PrintHits = 0); 
    
    static void printCaloProtoCluster(const mu2e::CaloProtoCluster* Cluster, const char* Opt = "");

    static void printCaloProtoClusterCollection(const mu2e::CaloProtoClusterCollection* Coll);
  };
}
#endif
