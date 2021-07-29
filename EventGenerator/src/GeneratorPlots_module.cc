// Module to plot generator output for any generator
// S. Middleton 2021

#include "Offline/MCDataProducts/inc/StageParticle.hh"
#include "Offline/MCDataProducts/inc/GenId.hh"
#include "Offline/DataProducts/inc/VirtualDetectorId.hh"

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
        fhicl::Atom<art::InputTag> genTag{Name("StageParticleCollection"), Comment("gen particle info")};
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
      const StageParticleCollection *_gencol;
      TTree* _Ntup;

      Int_t _genPdgId;
      Int_t _genCrCode;

      Float_t _genmomX;
      Float_t _genmomY;
      Float_t _genmomZ;
      Float_t _genStartX;
      Float_t _genStartY;
      Float_t _genStartZ;
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
    _Ntup->Branch("genMomX",      &_genmomX,      "genMomX/F");
    _Ntup->Branch("genMomY",      &_genmomY,      "genMomY/F");
    _Ntup->Branch("genMomZ",      &_genmomZ,      "genMomZ/F");
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
	
  for (unsigned int i=0; i <_gencol->size(); ++i)
  {
    StageParticle const& gen = (*_gencol)[i];
    _genPdgId   = gen.pdgId();
    _genCrCode  = gen.creationCode();
    _genmomX    = gen.momentum().x();
    _genmomY    = gen.momentum().y();
    _genmomZ    = gen.momentum().z();
    art::Ptr<mu2e::SimParticle> parent = gen.parent();
    _genStartX  = parent->endPosition().x();
    _genStartY  = parent->endPosition().y();
    _genStartZ  = parent->endPosition().z();
    _genStartT  = gen.time();
    _Ntup->Fill();
  } 
}


  bool GeneratorPlots::findData(const art::Event& evt){
    _gencol=0;
    auto genpart = evt.getValidHandle<StageParticleCollection>(_genTag);
    _gencol = genpart.product();
    return  _gencol!=0 ;
  }

 void GeneratorPlots::endJob(){}

}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::GeneratorPlots);

