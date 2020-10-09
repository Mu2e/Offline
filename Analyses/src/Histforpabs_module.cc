//
// A plugin to test using root interactively.
//
//
// Original author Rob Kutschke
//


// C++ includes.
#include <iostream>
#include <string>
#include <cmath>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

// Mu2e includes.
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"

// Root includes.
#include "TApplication.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TH1F.h"

using namespace std;

namespace mu2e {

  class Straw;

  class Histforpabs : public art::EDAnalyzer {
  public:

    explicit Histforpabs(fhicl::ParameterSet const& pset);
    virtual ~Histforpabs() { }

    virtual void beginJob();

    // This is called for each event.
    void analyze(const art::Event& e );

  private:

    // Start: run time parameters

    // Module label of the g4 module that made the hits.
    std::string _g4ModuleLabel;

    // Name of the tracker StepPoint collection
    std::string _trackerStepPoints;

    // Cut on the minimum energy.
    double _minimumEnergy;

    // Number of events analyzed.
    int _nAnalyzed;

    TH1F* _hEnergyat0;
    TH1F* _hEnergyat1;
    TH1F* _hEnergysim;

    void FillHistograms(const art::Event& event);
  };

  Histforpabs::Histforpabs(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),

    // Run time parameters
    _g4ModuleLabel(pset.get<string>("g4ModuleLabel")),
    _trackerStepPoints(pset.get<string>("trackerStepPoints","tracker")),
    _minimumEnergy(pset.get<double>("minimumEnergy")),

    // Histograms
    _hEnergyat0(0),
    _hEnergyat1(0),
    _hEnergysim(0){}

  void Histforpabs::beginJob( ){

    art::ServiceHandle<art::TFileService> tfs;
    _hEnergyat0 = tfs->make<TH1F>( "hEnergyat0", "Energy Deposited before 1st straw hist", 80, 102., 106.);
    _hEnergyat1 = tfs->make<TH1F>( "hEnergyat1", "Energy Deposited after 1st straw hist", 80, 102., 106.);
    _hEnergysim = tfs->make<TH1F>( "hEnergysim", "Sim particle energy", 80, 102.0, 106.0);

  }


  void Histforpabs::analyze(const art::Event& event ) {
    ++_nAnalyzed;
    FillHistograms(event);
  } // end analyze





  void Histforpabs::FillHistograms(const art::Event& event){

    art::Handle<StepPointMCCollection> hits;
    event.getByLabel(_g4ModuleLabel,_trackerStepPoints,hits);

    art::Handle<SimParticleCollection> simParticles;
    event.getByLabel(_g4ModuleLabel, simParticles);

    for( size_t i=0; i<hits->size(); ++i ){

      const StepPointMC& hit = (*hits)[i];
      // Skip hits with low pulse height.
      if ( hit.eDep() < _minimumEnergy ) continue;

      const CLHEP::Hep3Vector& mom = hit.momentum();
      SimParticleCollection::key_type trackId = hit.trackId();

      // Fill some histograms
      SimParticle const& simfm = simParticles->at(trackId);
      double trkrestm = simfm.endMomentum().e();
      if(i==0)_hEnergyat0->Fill(sqrt(mom.mag2() + trkrestm*trkrestm));
      if(i==1)_hEnergyat1->Fill(sqrt(mom.mag2() + trkrestm*trkrestm));

    } // end loop over hits.

    for ( SimParticleCollection::const_iterator i=simParticles->begin();
          i!=simParticles->end(); ++i ){

      SimParticle const& sim = i->second;
      if(!sim.madeInG4())_hEnergysim->Fill(sim.startMomentum().e());
    }
  } // end FillHistograms




}  // end namespace mu2e

using mu2e::Histforpabs;
DEFINE_ART_MODULE(Histforpabs);
