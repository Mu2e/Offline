// Module to plot generator output for any generator
// S. Middleton 2021

#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/GenId.hh"
#include "Offline/DataProducts/inc/XYZVec.hh"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Provenance.h"
#include "cetlib_except/exception.h"
#include "Offline/GeneralUtilities/inc/ParameterSetHelpers.hh"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <iostream>
#include <vector>

#include "TTree.h"
using namespace std;
namespace mu2e {

  class GeneratorPlots : public art::EDAnalyzer {
     
     public:

      struct Config {
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;
        fhicl::Atom<int> diagLevel{Name("diagLevel"),Comment("diag level"),0};
        fhicl::Atom<art::InputTag> genTag{Name("SimParticleCollection"), Comment("gen particle info")};
      };
      typedef art::EDAnalyzer::Table<Config> Parameters;

      explicit GeneratorPlots(const Parameters& conf);
      virtual ~GeneratorPlots() {}


      virtual void beginJob();
      virtual void endJob();
      virtual void analyze(const art::Event& e) override;

     private:
      Config _conf;
      int _diagLevel;
      art::InputTag _genTag;
      const SimParticleCollection *_gencol;
      TTree* _Ntup;

      Int_t _genPdgId;
      Int_t _genCrCode;

      Float_t _genmomX;
      Float_t _genmomY;
      Float_t _genmomZ;
      Float_t _genStartX;
      Float_t _genStartY;
      Float_t _genStartZ;
      XYZVec _genmom;
      XYZVec _genpos;
      Float_t _genStartT;
      bool findData(const art::Event& evt);
      void GetGenPartInfo(const art::Event& evt);
};


  GeneratorPlots::GeneratorPlots(const Parameters& conf):
  art::EDAnalyzer(conf),
  _diagLevel(conf().diagLevel()),
  _genTag(conf().genTag())
  {}

  void GeneratorPlots::beginJob(){
    art::ServiceHandle<art::TFileService> tfs;
    _Ntup  = tfs->make<TTree>("GenAna", "GenAna");
    _Ntup->Branch("genId",        &_genPdgId,     "genId/I");
    _Ntup->Branch("genCrCode",    &_genCrCode,    "genCrCode/I");
    _Ntup->Branch("genmom",       &_genmom,       "genmom/F");
    _Ntup->Branch("genMomX",      &_genmomX,      "genMomX/F");
    _Ntup->Branch("genMomY",      &_genmomY,      "genMomY/F");
    _Ntup->Branch("genMomZ",      &_genmomZ,      "genMomZ/F");
    _Ntup->Branch("genPos",       &_genpos,       "genpos/F");
    _Ntup->Branch("genStartX",    &_genStartX,    "genStartX/F");
    _Ntup->Branch("genStartY",    &_genStartY,    "genStartY/F");
    _Ntup->Branch("genStartZ",    &_genStartZ,    "genStartZ/F");
    _Ntup->Branch("genStartT",    &_genStartT,    "genStartT/F");
  }


  void GeneratorPlots::analyze(const art::Event& event) {

    if(!findData(event)) 
      throw cet::exception("RECO")<<"No data in  event"<< endl; 

    GetGenPartInfo(event);
}

void GeneratorPlots::GetGenPartInfo(const art::Event& evt){
	
  cet::map_vector<mu2e::SimParticle>::const_iterator iter;
  for(iter=_gencol->begin(); iter!=_gencol->end(); iter++)
  {
    const mu2e::SimParticle& particle = iter->second;
    _genPdgId   = particle.pdgId();
    _genCrCode  = particle.creationCode();
    _genmom = particle.startMomentum();
    _genpos = particle.startPosition();
    _genmomX    = particle.startMomentum().x();
    _genmomY    = particle.startMomentum().y();
    _genmomZ    = particle.startMomentum().z();
    _genStartX  = particle.startPosition().x();
    _genStartY  = particle.startPosition().y();
    _genStartZ  = particle.startPosition().z();
    _genStartT  = particle.startGlobalTime();
    _Ntup->Fill();
  } 
}


  bool GeneratorPlots::findData(const art::Event& evt){
    _gencol=0;
    auto genpart = evt.getValidHandle<SimParticleCollection>(_genTag);
    _gencol = genpart.product();
    return  _gencol!=0 ;
  }

 void GeneratorPlots::endJob(){}

}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::GeneratorPlots);

