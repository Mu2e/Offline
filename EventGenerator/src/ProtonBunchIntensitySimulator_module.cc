//
//  Module to simulate the microbunch-to-microbunch intensity fluctuations of the protons hitting the
//  primary target.  This module does not describe the time structure of the proton bunches, which is assumed
//  to be independent of the intensity.  The data product of this module is used by downstream simulation
//  modules, in particular, all generators which start from stopped muons  need to use the same (coherent)
//  intensity fluctuations for the rates of signal and background processes to be consistent.
//  Several possible models of intensity distributions can be selected
//  Original author: David Brown (LBNL) 19 May 2015
//
// art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
// general Mu2e includes
#include "SeedService/inc/SeedService.hh"
// CLHEP
#include "CLHEP/Random/RandFlat.h"
// data products produced by this module
#include "MCDataProducts/inc/ProtonBunchIntensity.hh"
#include <iostream>

namespace mu2e {
  class ProtonBunchIntensitySimulator : public art::EDProducer {
    public:
      enum ProtonIntensityModel {constant=0,flat};
      explicit ProtonBunchIntensitySimulator(const fhicl::ParameterSet& pset);
      virtual void produce(art::Event& event);
    private:
      ProtonIntensityModel _model;
      double _nmean; // mean number of protons/microbunch hitting the target
      double _width; // fractional full width of the flat distribution
      int _printLevel; // level of diagnostic printout
      double _flimits[2]; // cache of the limits for the flat distribution
      art::RandomNumberGenerator::base_engine_t& _engine;
      CLHEP::RandFlat _randflat; // flat generator
  };


  ProtonBunchIntensitySimulator::ProtonBunchIntensitySimulator(const fhicl::ParameterSet& pset) :
    _model(static_cast<ProtonIntensityModel>(pset.get<int>("IntensityModel",flat))),
    _nmean(pset.get<double>("MeanNumberOfProtonsPerMicrobunch")),
    _width(pset.get<double>("FullRelativeWidth",0.5)),
    _printLevel(pset.get<int>("PrintLevel",0)),
    _engine(createEngine( art::ServiceHandle<SeedService>()->getSeed())),
    _randflat(_engine)
  {
  // setup limits
    _flimits[0] = _nmean*(1.0 - 0.5*_width);
    _flimits[1] = _nmean*(1.0 + 0.5*_width);
  // setup random generator
    produces<mu2e::ProtonBunchIntensity>();
    if(_printLevel > 0){
      if(_model == constant)
	std::cout << "Generating proton bunches of constant intensity = " << _nmean << std::endl;
      else if(_model == flat)
	std::cout << "Generating proton bunches with flat intensity between " << _flimits[0] 
	<< " and " << _flimits[1] << std::endl;
      else
	std::cout <<"Error: unknown model " << _model << std::endl;
    }
  }

  void ProtonBunchIntensitySimulator::produce(art::Event& event) {
  // calculate intensity
    double fintensity;
    switch(_model) {
      case constant : default:
	fintensity = _nmean;
	break;
      case flat:
	fintensity = _randflat.fire(_flimits[0],_flimits[1]);
	break;
    }
    // convert to nearest ingeger
    unsigned intensity = static_cast<unsigned>(rint(fintensity));
    // create proton intensity object and put it into the event
    std::unique_ptr<mu2e::ProtonBunchIntensity> pbi ( new ProtonBunchIntensity(intensity, _nmean) );
    if(_printLevel > 1){
      std::cout << "Generating " << pbi->intensity() << " protons in this microbunch, from a mean of " 
      << pbi->meanIntensity() << std::endl;
    }
    event.put(std::move(pbi));
  }
}

DEFINE_ART_MODULE(mu2e::ProtonBunchIntensitySimulator);

