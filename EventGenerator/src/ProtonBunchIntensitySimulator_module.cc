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
#include "MCDataProducts/inc/EventWeight.hh"


namespace mu2e {
  class ProtonBunchIntensitySimulator : public art::EDProducer {
    public:
      enum ProtonIntensityModel {constant=0,flat};
      explicit ProtonBunchIntensitySimulator(const fhicl::ParameterSet& pset);
      virtual void produce(art::Event& event);
    private:
      ProtonIntensityModel _model;
      double _nmedian; // median number of protons/microbunch hitting the target
      double _width; // fractional full width of the flat distribution
      double _flimits[2]; // cache of the limits for the flat distribution
      art::RandomNumberGenerator::base_engine_t& _engine;
      CLHEP::RandFlat _randflat; // flat generator
  };


  ProtonBunchIntensitySimulator::ProtonBunchIntensitySimulator(const fhicl::ParameterSet& pset) :
    _model(static_cast<ProtonIntensityModel>(pset.get<int>("IntensityModel",flat))),
    _nmedian(pset.get<double>("MedianNumberOfProtonsPerMicrobunch")),
    _width(pset.get<double>("RelativeWidth",0.5)),
    _engine(createEngine( art::ServiceHandle<SeedService>()->getSeed())),
    _randflat(_engine)
  {
  // setup limits
    _flimits[0] = _nmedian + 0.5*_width;
    _flimits[1] = _nmedian - 0.5*_width;
  // setup random generator
    produces<mu2e::ProtonBunchIntensity>();
    produces<mu2e::EventWeight>();
  }

  void ProtonBunchIntensitySimulator::produce(art::Event& event) {
  // calculate intensity
    double fintensity;
    switch(_model) {
      case constant : default:
	fintensity = _nmedian;
	break;
      case flat:
	fintensity = _randflat.fire(_flimits[0],_flimits[1]);
	break;
    }
    // convert to nearest ingeger
    unsigned intensity = static_cast<unsigned>(rint(fintensity));
    // create proton intensity object and put it into the event
    std::unique_ptr<mu2e::ProtonBunchIntensity> pbi ( new ProtonBunchIntensity(intensity) );
    event.put(std::move(pbi));
    // downstream modules need to weight any process that depends on the proton bunch intensity, even
    // if they are only generating 1/event (like conversion electrons), as the probability of producing
    // that one event scales with proton intensity
    double weight = fintensity/_nmedian;
    std::unique_ptr<mu2e::EventWeight> evtwt ( new EventWeight(weight) );
    event.put(std::move(evtwt));
  }
}

DEFINE_ART_MODULE(mu2e::ProtonBunchIntensitySimulator);

