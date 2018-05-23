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
#include "art/Framework/Services/Optional/TFileService.h"
#include "cetlib_except/exception.h"
// root
#include "TH1F.h"
// general Mu2e includes
#include "SeedService/inc/SeedService.hh"
// CLHEP
#include "CLHEP/Random/RandFlat.h"
// data products produced by this module
#include "MCDataProducts/inc/ProtonBunchIntensity.hh"
#include <iostream>
#include <random>

namespace mu2e {
// wrapper around art engine that satisfies C++ URBG standard
  class artURBG {
    public:
      artURBG(art::RandomNumberGenerator::base_engine_t& engine) : _engine(engine) {}
      typedef unsigned int result_type;
      result_type min() { return 0; }
      result_type max() { return UINT_MAX; }
      result_type operator() () { return _engine.operator unsigned int(); }
    private:
      art::RandomNumberGenerator::base_engine_t& _engine; // CLHEP ENGINE
  };

  class ProtonBunchIntensitySimulator : public art::EDProducer {
    public:
      enum ProtonIntensityModel {constant=0,flat,lognorm};
      explicit ProtonBunchIntensitySimulator(const fhicl::ParameterSet& pset);
      virtual void beginJob();
      virtual void produce(art::Event& event);
    private:
      ProtonIntensityModel _model;
      double _nmean; // mean number of protons/microbunch hitting the target
      double _width; // fractional full width of the flat distribution
      double _sig, _mu, _lnmean; // lognormal parameters
      int _diag; // level of diag histograms
      int _printLevel; // level of diagnostic printout
      double _flimits[2]; // cache of the limits for the flat distribution
      art::RandomNumberGenerator::base_engine_t& _engine;
      CLHEP::RandFlat _randflat; // flat generator
      std::lognormal_distribution<double> _lognd;
      artURBG _urbg;
      // diagnostics
      TH1F *_pbi, *_pbir;
  };


  ProtonBunchIntensitySimulator::ProtonBunchIntensitySimulator(const fhicl::ParameterSet& pset) :
    _model(static_cast<ProtonIntensityModel>(pset.get<int>("IntensityModel",flat))),
    _nmean(pset.get<double>("MeanNumberOfProtonsPerMicrobunch")),
    _width(pset.get<double>("FullRelativeWidth",0.5)),
    _diag(pset.get<int>("DiagLevel",0)),
    _printLevel(pset.get<int>("PrintLevel",0)),
    _engine(createEngine( art::ServiceHandle<SeedService>()->getSeed())),
    _randflat(_engine),
    _lognd(pset.get<double>("Lognormal_mu",-0.1005),pset.get<double>("Lognormal_sigma",0.3814)),
    _urbg(_engine)
  {
  // setup limits
    _flimits[0] = _nmean*(1.0 - 0.5*_width);
    _flimits[1] = _nmean*(1.0 + 0.5*_width);
  // calculate lognormal mean
  _lnmean = exp(_mu + 0.5*_sig*_sig);
  // setup random generator
    produces<mu2e::ProtonBunchIntensity>();
    if(_printLevel > 0){
      if(_model == constant)
	std::cout << "Generating proton bunches of constant intensity = " << _nmean << std::endl;
      else if(_model == flat)
	std::cout << "Generating proton bunches with flat intensity between " << _flimits[0] 
	<< " and " << _flimits[1] << std::endl;
      else if(_model == lognorm)
	std::cout << "Generating proton bunches with lognormal intensity sigma =  " << _lognd.s()
	<< " mu: = " << _lognd.m() << " mean = " << _lnmean << std::endl;
      else
        throw cet::exception("SIM")<<"mu2e::ProtonBunchIntensitySimulator: unknown proton bunch intensity model " << _model << std::endl;
    }
  }

  void ProtonBunchIntensitySimulator::beginJob(){

    if(_diag > 0){
      art::ServiceHandle<art::TFileService> tfs;
      _pbi = tfs->make<TH1F>("pbi","Proton Bunch Intensity",100,0,1.2e8);
      _pbir = tfs->make<TH1F>("pbir","Relative Proton Bunch Intensity",100,0,3.0);

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
      case lognorm:
	fintensity = _nmean*_lognd(_urbg)/_lnmean;
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
    // optional diagnostics
    if(_diag > 0){
      _pbi->Fill(fintensity);
      _pbir->Fill(fintensity/_nmean);
    }
    event.put(std::move(pbi));
  }
}

DEFINE_ART_MODULE(mu2e::ProtonBunchIntensitySimulator);

