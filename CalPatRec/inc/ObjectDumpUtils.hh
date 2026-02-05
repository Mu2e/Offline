//
#ifndef CalPatRec_ObjectDumpUtils_hh
#define CalPatRec_ObjectDumpUtils_hh

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


#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/RecoDataProducts/inc/CaloProtoCluster.hh"

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

  public:
//-----------------------------------------------------------------------------
// methods
//-----------------------------------------------------------------------------

    static void SetFlagBgrHitsModuleLabel(const char* Label) { _FlagBgrHitsModuleLabel = Label; }

    static void printEventHeader(const art::Event* Event, const char* Message = "");

    static void printCaloProtoCluster(const mu2e::CaloProtoCluster* Cluster, const char* Opt = "");

    static void printCaloProtoClusterCollection(const mu2e::CaloProtoClusterCollection* Coll);
  };
}
#endif
