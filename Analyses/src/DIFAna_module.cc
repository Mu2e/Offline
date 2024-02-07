//
// TODO:  This file will analyze the outgoing photons energy, momentum and position by looking at the indicent electron/positron
// Original author S. Middleton
//

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"

#include "Offline/ProditionsService/inc/ProditionsHandle.hh"

#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"


#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include <iostream>
#include <fstream>
#include <vector>

//ROOT
#include "TTree.h"

using namespace std;
namespace mu2e{

  class DIFAna : public art::EDAnalyzer {
	  public:
		  struct Config {
		  using Name=fhicl::Name;
		  using Comment=fhicl::Comment;
		  fhicl::Atom<int> diag{Name("diag"), Comment("Create diag histograms"),0};
	    fhicl::Atom<art::InputTag> SimToken{Name("SimParticleCollection"),Comment("tag for Sim collection")};
  };
	  typedef art::EDAnalyzer::Table<Config> Parameters;
	  explicit DIFAna(const Parameters& conf);
	  virtual void beginJob() override;
	  virtual void analyze(const art::Event& e);

  private:
	  Config _conf;
	  int _diagLevel;
	  art::InputTag _SimToken;
	  const SimParticleCollection* _SimCol;
	  
    TTree* genTree;
    //Float_t _maxr;
    Float_t _momT;
    Float_t _posT;
    Float_t _cosTheta;
    Float_t _time;

	  };

  DIFAna::DIFAna(const Parameters& conf) :
  art::EDAnalyzer(conf),
  _diagLevel(conf().diag()),
  _SimToken(conf().SimToken())
  {}

  void DIFAna::beginJob() { 
	  art::ServiceHandle<art::TFileService> tfs;
      genTree  = tfs->make<TTree>("GenAna", "GenAna");
      //genTree->Branch("maxr", &_maxr, "maxr/F");   
      genTree->Branch("momT", &_momT, "momT/F"); 
      genTree->Branch("posT", &_posT, "posT/F");
      genTree->Branch("cosTheta", &_cosTheta, "cosTheta/F");
      genTree->Branch("time", &_time, "time/F");
  }
  void DIFAna::analyze(const art::Event& event) {

    //------------SimParticles-------------//
    auto sH = event.getValidHandle<mu2e::SimParticleCollection>(_SimToken);
    _SimCol = sH.product();

    for ( SimParticleCollection::const_iterator i=_SimCol->begin(); i!=_SimCol->end(); ++i ){
      SimParticle const& sim = i->second;
      _momT = sim.startMomentum().rho();
      _posT = sim.startPosition().rho();
      _cosTheta = cos(atan2(_momT,sim.startMomentum().z()));
      _time = sim.time();
      genTree->Fill();
      }

  }
}//end mu2e namespace
using mu2e::DIFAna;
DEFINE_ART_MODULE(DIFAna)
