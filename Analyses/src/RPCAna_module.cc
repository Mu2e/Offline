//Author: S Middleton
//Date April 2020
//Purpose: Basic Skeleton for GenParticle Aanlyzer
#include "CLHEP/Units/SystemOfUnits.h"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GlobalConstantsService/inc/unknownPDGIdName.hh"

#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
//#include "MCDataProducts/inc/StatusG4.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Provenance.h"
#include "cetlib_except/exception.h"
#include "GeneralUtilities/inc/ParameterSetHelpers.hh"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "TNtuple.h"

#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <memory>
#include <vector>

using namespace std;
using CLHEP::Hep3Vector;
using CLHEP::keV;

namespace mu2e {

  class RPCAna : public art::EDAnalyzer {

     public:

     struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<int> diagLevel{Name("diagLevel"),Comment("diag level"),0};
      fhicl::Atom<int> mcdiag{Name("mcdiag"),Comment("mc diag level"),0};
      fhicl::Atom<art::InputTag> genTag{Name("GenParticleCollection"), Comment("gen particle info")};
      fhicl::Atom<art::InputTag> stepTag{Name("StepPointMCCollection"), Comment("step point info")};
      fhicl::Atom<art::InputTag> simTag{Name("SimParticleCollection"), Comment("sim particle info")};
      };
       typedef art::EDAnalyzer::Table<Config> Parameters;

       explicit RPCAna(const Parameters& conf);
       virtual ~RPCAna() {}


       virtual void beginJob();
       virtual void endJob();
	     virtual void analyze(const art::Event& e) override;

     private:
       Config _conf;
       int _diagLevel;
       int _nProcess;
       int _mcdiag;

       art::InputTag _genTag;
       art::InputTag _stepTag;
       art::InputTag _simTag;
      
       const GenParticleCollection *_gencol;
       const SimParticleCollection *_simcol;
       const StepPointMCCollection *_stepcol;

       TTree* _Ntup;
       Int_t _evt, _run, _nSim, _nGen,_genPdgId, _genCrCode;
       Float_t _genmomX,_genmomY,_genmomZ, _genmomTot, _genStartX,_genStartY,_genStartZ,_genStartT;
       
       bool findData(const art::Event& evt);

};


  RPCAna::RPCAna(const Parameters& conf):
    art::EDAnalyzer(conf),
    _diagLevel(conf().diagLevel()),
    _mcdiag(conf().mcdiag()),
    _genTag(conf().genTag()),
    _stepTag(conf().stepTag()),
    _simTag(conf().simTag())
{}

  void RPCAna::beginJob(){

    art::ServiceHandle<art::TFileService> tfs;
    _Ntup  = tfs->make<TTree>("RPCAna", "RPCAna");
    _Ntup->Branch("evt",          &_evt ,        "evt/I");
    _Ntup->Branch("run",          &_run ,        "run/I");
    _Ntup->Branch("nGen",         &_nGen ,        "nGen/I");
    _Ntup->Branch("genId",        &_genPdgId,     "genId/I");
    _Ntup->Branch("genCrCode",    &_genCrCode,    "genCrCode/I");
    _Ntup->Branch("genMomX",      &_genmomX,      "genMomX/F");
    _Ntup->Branch("genMomY",      &_genmomY,      "genMomY/F");
    _Ntup->Branch("genMomZ",      &_genmomZ,      "genMomZ/F");
    _Ntup->Branch("genMomTot",    &_genmomTot,      "genMomTot/F");
    _Ntup->Branch("genStartX",    &_genStartX,    "genStartX/F");
    _Ntup->Branch("genStartY",    &_genStartY,    "genStartY/F");
    _Ntup->Branch("genStartZ",    &_genStartZ,    "genStartZ/F");
    _Ntup->Branch("genStartT",    &_genStartT,    "genStartT/F");

  }


  void RPCAna::analyze(const art::Event& event) {

   _evt = event.id().event();
   _run = event.run();

   if(!findData(event))
  		throw cet::exception("RECO")<<"No data in  event"<< endl;

  ++_nProcess;
  if (_nProcess%1000==0) std::cout<<"Processing event "<<_nProcess<<std::endl;
  _nGen = _gencol->size();
  for (unsigned int i=0; i <_gencol->size(); ++i)
  {
     GenParticle const& gen = (*_gencol)[i];
    _genPdgId   = gen.pdgId();
    _genCrCode  = gen.generatorId().id();
    _genmomX    = gen.momentum().vect().x();
    _genmomY    = gen.momentum().vect().y();
    _genmomZ    = gen.momentum().vect().z();
    _genmomTot  = gen.momentum().vect().x()+gen.momentum().vect().y()+ gen.momentum().vect().z();
    _genStartX  = gen.position().x();
    _genStartY  = gen.position().y();
    _genStartZ  = gen.position().z();
    _genStartT  = gen.time();
  }
   _Ntup->Fill();
}

bool RPCAna::findData(const art::Event& evt){

  _gencol=0;
  auto genpart = evt.getValidHandle<GenParticleCollection>(_genTag);
  _gencol = genpart.product();

  _simcol=0;
  auto simpart = evt.getValidHandle<SimParticleCollection>(_simTag);
  _simcol = simpart.product();

  _stepcol=0;
  auto step = evt.getValidHandle<StepPointMCCollection>(_stepTag);
  _stepcol = step.product();

	return  _gencol!=0 and _stepcol!=0 and _simcol!=0; //Add in other collections
}

 void RPCAna::endJob(){}

}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::RPCAna);
