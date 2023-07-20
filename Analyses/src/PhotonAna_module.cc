//
// TODO:  This file will analyze the outgoing photons energy, momentum and position by looking at the indicent electron/positron
// Original author S. Middleton and H. Jafree
//

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"

#include "Offline/ProditionsService/inc/ProditionsHandle.hh"

#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"

#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <utility>
#include <functional>
#include <float.h>
#include <vector>
#include <map>

//ROOT
#include "TStyle.h"
#include "Rtypes.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TTree.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TProfile.h"
using namespace std;
namespace mu2e{

  class PhotonAna : public art::EDAnalyzer {
    public:
      struct Config {
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;
        fhicl::Atom<int> diag{Name("diag"), Comment("Create diag histograms"),0};
        fhicl::Atom<art::InputTag> KalToken{Name("KalSeedCollection"),Comment("tag for kal seed collection")};
      };
      typedef art::EDAnalyzer::Table<Config> Parameters;
        explicit PhotonAna(const Parameters& conf);
      virtual ~PhotonAna(){};
      virtual void beginJob() override;
      virtual void analyze(const art::Event& e);
      
    private:

      Config _conf;
      int _diag;
      art::InputTag _KalToken;
      const KalSeedCollection* _KalCol;

      TTree* _photon_analyzer; 
      Float_t _pathlength;

  };

  PhotonAna::PhotonAna(const Parameters& conf) :
  art::EDAnalyzer(conf),
  _diag(conf().diag()),
  _KalToken(conf().KalToken())
  {}
  

  void PhotonAna::beginJob() { //TODO - can add TTree and THistograms here if required
      if(_diag > 0){
      art::ServiceHandle<art::TFileService> tfs;
      _photon_analyzer=tfs->make<TTree>("photon_analyzer"," Diagnostics for Photon Conversion Track Fitting"); 
      _photon_analyzer->Branch("PathLength", &_pathlength, "PathLength/F"); 
    }
  }

  using LHPT = KinKal::PiecewiseTrajectory<KinKal::LoopHelix>;
  void PhotonAna::analyze(const art::Event& event) {
    auto kalH = event.getValidHandle<KalSeedCollection>(_KalToken);
    _KalCol = kalH.product();
    for(unsigned int k = 0; k < _KalCol->size(); k++){
      KalSeed kseed = (*_KalCol)[k];
      if(kseed.loopHelixFit()){
        std::unique_ptr<LHPT> trajectory = kseed.loopHelixFitTrajectory();
        double t1 = trajectory->range().begin();
        double x1 = trajectory->position3(t1).x();
        double y1 = trajectory->position3(t1).y();
        double z1 = trajectory->position3(t1).z();
        _pathlength = sqrt((x1)*(x1)+(y1)*(y1)+(z1)*(z1));
        _photon_analyzer->Fill(); 
        }
      }
    } 

}//end mu2e namespace
using mu2e::PhotonAna;
DEFINE_ART_MODULE(PhotonAna);
