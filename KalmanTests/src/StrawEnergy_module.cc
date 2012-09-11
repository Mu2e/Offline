//
//  Study energy desposited in straws 
//
// $Id: StrawEnergy_module.cc,v 1.3 2012/09/11 21:46:41 brownd Exp $
// $Author: 
// $Date: 2012/09/11 21:46:41 $
//
// Original author David Brown
//
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "ConditionsService/inc/GlobalConstantsHandle.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "MCDataProducts/inc/VirtualDetectorId.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "KalmanTests/inc/KalFitMC.hh"
#include "TTree.h"
#include "TBranch.h"
#include <cmath>
#include <iostream>
#include <string>

using namespace std;

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
    std::string _makerModuleLabel;
    TTree* _estraw;
// branch variables
    Int_t _device, _sector, _layer, _straw;
    Int_t _mcpdg, _mcgen, _mcproc;
    Int_t _mcppdg;
    threevec _mcpspos, _mcpsmom;
    Float_t _mcpke, _mcptime;
    Float_t _edep, _mce;
    Float_t _hittime;
    Float_t _strawx, _strawrho, _strawz;
  // helper functions
    art::Ptr<SimParticle> findPrimary(PtrStepPointMCVector const& mcptr) const;
    art::Ptr<SimParticle> findPrimary(StepPointMC const& mcstep) const;
};


  StrawEnergy::StrawEnergy(fhicl::ParameterSet const& pset) :
    // Run time parameters
    _diagLevel(pset.get<int>("diagLevel",0)),
    _makerModuleLabel(pset.get<string>("makerModuleLabel","makeSH"))
  {

  }

  void StrawEnergy::beginJob(){
    // Get access to the TFile service.
    art::ServiceHandle<art::TFileService> tfs;
    _estraw=tfs->make<TTree>("estraw","straw energy");
    _estraw->Branch("device",&_device,"device/I");
    _estraw->Branch("sector",&_sector,"sector/I");
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
    const Tracker& tracker = getTrackerOrThrow();

    art::Handle<StrawHitCollection> pdataHandle;
    event.getByLabel(_makerModuleLabel,pdataHandle);
    StrawHitCollection const* hits = pdataHandle.product();

    art::Handle<PtrStepPointMCVectorCollection> mcptrHandle;
    event.getByLabel(_makerModuleLabel,"StrawHitMCPtr",mcptrHandle);
    PtrStepPointMCVectorCollection const* hits_mcptr = mcptrHandle.product();
    // Loop over all hits
    for ( size_t ihit=0; ihit<hits->size(); ++ihit ){

      StrawHit const& hit(hits->at(ihit));
      _edep = hit.energyDep();
      _hittime = hit.time();
      // Get the hit poosition information.  This requires MC truth
      CLHEP::Hep3Vector mcpos;
// there can be more than 1 StepPointMC for each straw hit
      PtrStepPointMCVector const& mcptr(hits_mcptr->at(ihit));
// summarize direct MC truth information
      std::vector<TrkSum> mcsum;
      KalFitMC::fillMCHitSum(mcptr,mcsum);
      _mcpdg = mcsum[0]._pdgid;
      _mcgen = mcsum[0]._gid;
      _mcproc = mcsum[0]._pid;
      _mce = mcsum[0]._esum;
// follow this StepPointMC SimParticle back to the 'primary', and record it's info
      art::Ptr<SimParticle> primary = findPrimary(mcptr);
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
      const Straw& straw = tracker.getStraw( hit.strawIndex() );
      _device = straw.id().getDevice();
      _sector = straw.id().getSector();
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
// find primary particle from a set of StepPointMCs
  art::Ptr<SimParticle>
  StrawEnergy::findPrimary(PtrStepPointMCVector const& mcptr) const {
    Hep3Vector pos;
    double energy(0.0);
    std::map< art::Ptr<SimParticle>, double> spmap;
    for (size_t imc = 0; imc < mcptr.size(); ++imc) {
      StepPointMC const& mchit = *mcptr[imc];
      pos += mchit.position()*mchit.eDep();
      energy += mchit.eDep();
      art::Ptr<SimParticle> const& primary = findPrimary(mchit);
      if(primary.isNonnull()){
	std::map< art::Ptr<SimParticle>, double>::iterator ifnd = spmap.find(primary);
	if(ifnd != spmap.end()){
	  // already found this primary: add this energy
	  ifnd->second += mchit.eDep();
	} else {
	  // new primary: initialize with this energy
	  spmap[primary] = mchit.eDep();
	}
      }
    }
    if(energy >0.0)pos *= 1.0/energy;
    // primary particle information
    std::map< art::Ptr<SimParticle>, double>::iterator ibest = spmap.begin();
    for(std::map< art::Ptr<SimParticle>, double>::iterator isp = spmap.begin()++; isp != spmap.end(); ++isp){
      if(isp->second > ibest->second)ibest = isp;
    }
    if(ibest != spmap.end())
      return ibest->first;
    else
      return art::Ptr<SimParticle>();
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
