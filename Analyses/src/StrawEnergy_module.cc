//
//  Study energy desposited in straws 
//
// $Id: StrawEnergy_module.cc,v 1.6 2014/08/22 19:55:50 brownd Exp $
// $Author: 
// $Date: 2014/08/22 19:55:50 $
//
// Original author David Brown
//
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
#include "DataProducts/inc/XYZVec.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "TrkDiag/inc/KalDiag.hh"
#include "TTree.h"
#include "TBranch.h"
#include <cmath>
#include <iostream>
#include <string>

using namespace std;
using CLHEP::Hep3Vector;

namespace mu2e {
 // simple structs

  class StrawEnergy : public art::EDAnalyzer {
  public:

    explicit StrawEnergy(fhicl::ParameterSet const& pset);
    virtual ~StrawEnergy() { }

    virtual void beginJob();
    virtual void endJob();

    // This is called for each event.
    virtual void analyze(const art::Event& e);

  private:
    int _diagLevel;
    std::string _makerModuleLabel,_mcdigislabel;
    TTree* _estraw;
// branch variables
    Int_t _plane, _panel, _layer, _straw;
    Int_t _mcpdg, _mcgen, _mcproc;
    Int_t _mcppdg;
    XYZVec _mcpspos, _mcpsmom;
    Float_t _mcpke, _mcptime;
    Float_t _edep, _mce;
    Float_t _hittime;
    Float_t _strawx, _strawrho, _strawz;
  // helper functions
    art::Ptr<SimParticle> findPrimary(StepPointMC const& mcstep) const;
};


  StrawEnergy::StrawEnergy(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    // Run time parameters
    _diagLevel(pset.get<int>("diagLevel",0)),
    _makerModuleLabel(pset.get<string>("makerModuleLabel","makeSH")),
    _mcdigislabel(pset.get<std::string>("StrawHitMCLabel","makeSH"))
  {

  }

  void StrawEnergy::beginJob(){
    // Get access to the TFile service.
    art::ServiceHandle<art::TFileService> tfs;
    _estraw=tfs->make<TTree>("estraw","straw energy");
    _estraw->Branch("plane",&_plane,"plane/I");
    _estraw->Branch("panel",&_panel,"panel/I");
    _estraw->Branch("layer",&_layer,"layer/I");
    _estraw->Branch("straw",&_straw,"straw/I");
    _estraw->Branch("mcpdg",&_mcpdg,"mcpdg/I");
    _estraw->Branch("mcgen",&_mcgen,"mcgen/I");
    _estraw->Branch("mcproc",&_mcproc,"mcproc/I");
    _estraw->Branch("mce",&_mce,"mce/F");
    _estraw->Branch("mcppdg",&_mcppdg,"mcppdg/I");
    _estraw->Branch("mcpspos",&_mcpspos,"pspx/F:pspy/F:pspz/F");
    _estraw->Branch("mcpsmom",&_mcpsmom,"psmx/F:psmy/F:psmz/F");
    _estraw->Branch("mcpke",&_mcpke,"mcpke/F");
    _estraw->Branch("mcptime",&_mcptime,"mcptime/F");
    _estraw->Branch("edep",&_edep,"edep/F");
    _estraw->Branch("hittime",&_hittime,"hittime/F");
    _estraw->Branch("strawx",&_strawx,"strawx/F");
    _estraw->Branch("strawrho",&_strawrho,"strawrho/F");
    _estraw->Branch("strawz",&_strawz,"strawz/F");
  }

  void StrawEnergy::endJob(){
  }

  void StrawEnergy::analyze(const art::Event& event) {
    GlobalConstantsHandle<ParticleDataTable> pdt;
    const Tracker& tracker = *GeomHandle<Tracker>();

    art::Handle<StrawHitCollection> pdataHandle;
    event.getByLabel(_makerModuleLabel,pdataHandle);
    StrawHitCollection const* hits = pdataHandle.product();

    art::Handle<StrawDigiMCCollection> mcdigisHandle;
    event.getByLabel(_mcdigislabel,"StrawHitMC",mcdigisHandle);
    StrawDigiMCCollection const* mcdigis = mcdigisHandle.product();

    // Loop over all hits
    for ( size_t ihit=0; ihit<hits->size(); ++ihit ){

      StrawHit const& hit(hits->at(ihit));
      _edep = hit.energyDep();
      _hittime = hit.time();
      // Get the hit poosition information.  This requires MC truth
      CLHEP::Hep3Vector mcpos;
// find MC truth information
      StrawDigiMC const& mcdigi = mcdigis->at(ihit);
      // use TDC channel 0 to define the MC match
      StrawEnd itdc(StrawEnd::cal);
      art::Ptr<StepPointMC> const& spmcp = mcdigi.stepPointMC(itdc);
      art::Ptr<SimParticle> const& spp = spmcp->simParticle();
      _mcpdg = spp->pdgId();
      _mcgen = -1;
      if(spp->genParticle().isNonnull())
	_mcgen = spp->genParticle()->generatorId().id();
      _mcproc = spp->creationCode();
      _mce = mcdigi.energySum();
// follow this StepPointMC SimParticle back to the 'primary', and record it's info
      art::Ptr<SimParticle> primary = findPrimary(*spmcp);
     if(primary.isNonnull()){
	_mcppdg = primary->pdgId();
	_mcpsmom = primary->startMomentum();
	_mcpspos = primary->startPosition();
	mcpos = primary->startPosition();
	double mass = pdt->particle(primary->pdgId()).ref().mass();
	_mcpke = primary->startMomentum().e()-mass;
	_mcptime = primary->startGlobalTime();
      } else {
	_mcppdg = 0;
	_mcpsmom = Hep3Vector(0.0,0.0,0.0);
	_mcpspos = Hep3Vector(0.0,0.0,0.0);
      	mcpos = Hep3Vector(0.0,0.0,0.0);
	_mcpke = 0.0;
	_mcptime = 0.0;
      }
// Get the straw information:
      const Straw& straw = tracker.getStraw( hit.strawId() );
      _plane = straw.id().getPlane();
      _panel = straw.id().getPanel();
      _layer = straw.id().getLayer();
      _straw = straw.id().getStraw();
      const CLHEP::Hep3Vector& mid   = straw.getMidPoint();
      const CLHEP::Hep3Vector& w     = straw.getDirection();
// define position using straw
      _strawx = w.dot(mcpos-mid);
      _strawrho = mid.perp();
      _strawz = mid.z();
// fill the tree
      _estraw->Fill();
    }
  }

// find the virtual detector hits and SimParticles for this event
  art::Ptr<SimParticle>
  StrawEnergy::findPrimary(StepPointMC const& mcstep) const {
  // start by finding the immediate parent
    art::Ptr<SimParticle> simp = mcstep.simParticle();
  // loop backwards until this is a primary
    while(simp.isNonnull() && !simp->isPrimary()){
      simp = simp->parent();
    }
    return simp;
  }
}

DEFINE_ART_MODULE(mu2e::StrawEnergy);
