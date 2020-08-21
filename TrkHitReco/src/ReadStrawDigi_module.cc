//
// This module transforms StrawDigi objects into StrawHit objects
//
//
// Original author David Brown, LBNL
// Merged with flag and position creation B. Echenard, CalTech
//
// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "RecoDataProducts/inc/StrawDigi.hh"
#include <iostream>

namespace mu2e {
  using namespace TrkTypes;

  class ReadStrawDigi : public art::EDProducer 
  {
    public:
      explicit ReadStrawDigi(fhicl::ParameterSet const& pset);
      virtual ~ReadStrawDigi(); 
      virtual void produce( art::Event& e);
    private:
      art::InputTag _sdtag;
  };

  ReadStrawDigi::ReadStrawDigi(fhicl::ParameterSet const& pset) :
      art::EDProducer{pset},
      _sdtag  (pset.get<art::InputTag>("StrawDigiCollection","makeSD"))
  { }

  ReadStrawDigi::~ReadStrawDigi() {}

  void ReadStrawDigi::produce(art::Event& event)
  {        
    auto sdH = event.getValidHandle<StrawDigiCollection>(_sdtag);
    const StrawDigiCollection& sdcol(*sdH);
     
    for (size_t isd=0;isd<sdcol.size();++isd) {
      const StrawDigi& digi = sdcol[isd];
      // access a few fields to force a read
      if(digi.strawId().asUint16() == 0xffff ||
	digi.TDC()[0] == 0xffff ||
	digi.TOT()[0] == 0xffff ||
	digi.adcWaveform()[0] == 0xffff)
	std::cout << "error " << std::endl;
    }
  }
}

using mu2e::ReadStrawDigi;
DEFINE_ART_MODULE(ReadStrawDigi);

