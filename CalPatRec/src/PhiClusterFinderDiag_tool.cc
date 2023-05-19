
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "fhiclcpp/ParameterSet.h"

#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "Offline/CalPatRec/inc/PhiClusterFinder_types.hh"
#include "Offline/Mu2eUtilities/inc/McUtilsToolBase.hh"
#include "Offline/Mu2eUtilities/inc/ModuleHistToolBase.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"

#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art/Framework/Principal/Event.h"

#include "TH1.h"

namespace mu2e {

  using namespace PhiClusterFinderTypes;
  using PhiClusterFinderTypes::Config;

  class PhiClusterFinderDiag : public mu2e::ModuleHistToolBase {

  public:

    enum {
      kNEventHistsSets = 1,
      kNPhiClusterHistsSets = 1
    };

    struct Hists {
       TH1F* nTimeClusters;
       TH1F* nTimeClustersNew;
       TH1F* nHits;
    };

  protected:
    Hists                      _hist;
    Data_t*                    _data;
    int                        _event_number;
    int                        _mcTruth;

  public:

    explicit PhiClusterFinderDiag(const Config& config);
    explicit PhiClusterFinderDiag(const fhicl::ParameterSet& PSet);
    ~PhiClusterFinderDiag();

  private:

    virtual int bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) override ;
    virtual int fillHistograms(void* Data, int Mode = -1) override ;
  };

//-----------------------------------------------------------------------------
  PhiClusterFinderDiag::PhiClusterFinderDiag(const Config& config) {
    printf(" PhiClusterFinderDiag::PhiClusterFinderDiag On\n");

  }

//-----------------------------------------------------------------------------
  PhiClusterFinderDiag::PhiClusterFinderDiag(const fhicl::ParameterSet& PSet) {
    printf(" PhiClusterFinderDiag::PhiClusterFinderDiag On\n");
  }

//-----------------------------------------------------------------------------
  PhiClusterFinderDiag::~PhiClusterFinderDiag() {}

//-----------------------------------------------------------------------------
  int PhiClusterFinderDiag::bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) {

    _hist.nTimeClusters    = Tfs->make<TH1F>("ntc"  , "number of time clusters"                      , 50, 0, 20);
    _hist.nTimeClustersNew   = Tfs->make<TH1F>("ntcnew", "number of phi clusters"             , 50, 0, 20);
    _hist.nHits       = Tfs->make<TH1F>("nhits"  , "number of hits"             , 100, 0, 100);
    return 0;

  }

//-----------------------------------------------------------------------------
// Mode is not used here
//-----------------------------------------------------------------------------
  int PhiClusterFinderDiag::fillHistograms(void* Data, int Mode) {
    _data = (Data_t*) Data;

    _hist.nTimeClusters ->Fill((float)_data->_tccol->size());
    _hist.nTimeClustersNew ->Fill((float)_data->_tccolnew->size());

    for (int i=0; i<(int)_data->_tccolnew->size(); i++) {
      float numHits = (float)_data->_tccolnew->at(i)._strawHitIdxs.size();
      // printf("Num hits = %f \n",numHits);
      _hist.nHits->Fill(numHits);
    }

    return 0;
  }

  DEFINE_ART_CLASS_TOOL(PhiClusterFinderDiag)

}
