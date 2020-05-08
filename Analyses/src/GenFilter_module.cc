
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "MCDataProducts/inc/CaloHitSimPartMCCollection.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"

// Framework includes.
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
// Root includes
//#include "TNtuple.h"
#include "TTree.h"

// Other includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ includes
#include <iostream>
#include <vector>

using namespace std;

namespace mu2e {

  class GenFilter: public art::EDFilter {

  public:
	struct Config {
		     using Name=fhicl::Name;
		     using Comment=fhicl::Comment;
	 	      fhicl::Atom<art::InputTag> genTag{Name("GenParticleCollection"), Comment("gen particle info")};
		     
	    	};
  
    typedef art::EDFilter::Table<Config> Parameters;
		explicit GenFilter(const Parameters& conf);
    virtual bool filter  (art::Event& event) override;
		virtual ~GenFilter() {}
   
  private:
    Config _conf;
  
    TTree* _Ntup;
    art::InputTag _genTag;
    Int_t _evt, _run, _nSim, _nGen,_genPdgId, _genCrCode;
    Float_t _genmomX,_genmomY,_genmomZ, _genmomTot, _genStartX,
    _genStartY,_genStartZ,_genStartT, _genKE;

    const GenParticleCollection *_gencol;
    bool findData(const art::Event& evt);
  };

  GenFilter::GenFilter(const Parameters& conf):
    art::EDFilter{conf},
    _genTag(conf().genTag())
{

    art::ServiceHandle<art::TFileService> tfs;
    
    _Ntup  = tfs->make<TTree>("ExampleTree", "ExampleTree");
    _Ntup->Branch("evt",          &_evt ,        "evt/I");
    _Ntup->Branch("run",          &_run ,        "run/I");
    _Ntup->Branch("nGen",         &_nGen ,        "nGen/I");
    _Ntup->Branch("genId",        &_genPdgId,     "genId/I");
    _Ntup->Branch("genCrCode",    &_genCrCode,    "genCrCode/I");
    _Ntup->Branch("genMomX",      &_genmomX,      "genMomX/F");
    _Ntup->Branch("genMomY",      &_genmomY,      "genMomY/F");
    _Ntup->Branch("genMomZ",      &_genmomZ,      "genMomZ/F");
    _Ntup->Branch("genMomTot",    &_genmomTot,    "genMomTot/F");
    _Ntup->Branch("genStartX",    &_genStartX,    "genStartX/F");
    _Ntup->Branch("genStartY",    &_genStartY,    "genStartY/F");
    _Ntup->Branch("genStartZ",    &_genStartZ,    "genStartZ/F");
    _Ntup->Branch("genStartT",    &_genStartT,    "genStartT/F");
    _Ntup->Branch("genKE",        &_genKE,    "genKE/F");

}

  bool GenFilter::filter(art::Event& event) {
   
   if(!findData(event))
  		throw cet::exception("RECO")<<"No data in  event"<< endl;
  bool returnOK = true;
  for (unsigned int i=0; i <_gencol->size(); ++i)
  {
    GenParticle const& gen = (*_gencol)[i];
    
    
    if (!isfinite(KE)){ 
      _genPdgId   = gen.pdgId();
      _genCrCode  = gen.generatorId().id();
      _genmomX    = gen.momentum().vect().x();
      _genmomY    = gen.momentum().vect().y();
      _genmomZ    = gen.momentum().vect().z();
      _genmomTot  = sqrt(_genmomX*_genmomX+_genmomY*_genmomY+_genmomZ*_genmomZ);
      _genStartX  = gen.position().x();
      _genStartY  = gen.position().y();
      _genStartZ  = gen.position().z();
      _genStartT  = gen.time();
      _Ntup->Fill();
      returnOK = false;
    }
  }
   return returnOK;
 }

bool GenFilter::findData(const art::Event& evt){ 
  _gencol=0;
  auto genpart = evt.getValidHandle<GenParticleCollection>(_genTag);
  _gencol = genpart.product();
	return  _gencol!=0 ;
 }
}

using mu2e::GenFilter;
DEFINE_ART_MODULE(GenFilter);

