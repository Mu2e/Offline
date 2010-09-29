//
// A plugin to test using root interactively.
//
// $Id: Histforpabs_plugin.cc,v 1.2 2010/09/29 22:34:44 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/09/29 22:34:44 $
//
// Original author Rob Kutschke
//


// C++ includes.
#include <iostream>
#include <string>
#include <cmath>

// Framework includes.
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Services/interface/TFileService.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

// Mu2e includes.
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "ToyDP/inc/SimParticleCollection.hh"

// Root includes.
#include "TApplication.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TH1F.h"

using namespace std;

namespace mu2e {

  class Straw;

  class Histforpabs : public edm::EDAnalyzer {
  public:
    
    explicit Histforpabs(edm::ParameterSet const& pset);
    virtual ~Histforpabs() { }

    virtual void beginJob(edm::EventSetup const&);
 
    // This is called for each event.
    void analyze(const edm::Event& e, edm::EventSetup const&);

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

    void Histforpabs::FillHistograms(const edm::Event& event);
  };

  Histforpabs::Histforpabs(edm::ParameterSet const& pset) : 

    // Run time parameters
    _g4ModuleLabel(pset.getParameter<string>("g4ModuleLabel")),
    _trackerStepPoints(pset.getUntrackedParameter<string>("trackerStepPoints","tracker")),
    _minimumEnergy(pset.getParameter<double>("minimumEnergy")),
    
    // Histograms
    _hEnergyat0(0),
    _hEnergyat1(0),
    _hEnergysim(0){}

  void Histforpabs::beginJob(edm::EventSetup const& ){
    
    edm::Service<edm::TFileService> tfs;
    _hEnergyat0 = tfs->make<TH1F>( "hEnergyat0", "Energy Deposited before 1st straw hist", 80, 102., 106.);
    _hEnergyat1 = tfs->make<TH1F>( "hEnergyat1", "Energy Deposited after 1st straw hist", 80, 102., 106.);
    _hEnergysim = tfs->make<TH1F>( "hEnergysim", "Sim particle energy", 80, 102.0, 106.0);
    
  }


  void Histforpabs::analyze(const edm::Event& event, edm::EventSetup const&) {
    ++_nAnalyzed;
    FillHistograms(event);   
  } // end analyze





  void Histforpabs::FillHistograms(const edm::Event& event){

    edm::Handle<StepPointMCCollection> hits;
    event.getByLabel(_g4ModuleLabel,_trackerStepPoints,hits);

    edm::Handle<SimParticleCollection> simParticles;
    event.getByType(simParticles);

    for( size_t i=0; i<hits->size(); ++i ){

      const StepPointMC& hit = (*hits)[i];
      // Skip hits with low pulse height.
      if ( hit.eDep() < _minimumEnergy ) continue;
     
      const CLHEP::Hep3Vector& mom = hit.momentum();  int trackId = hit.trackId();

      // Fill some histograms
      SimParticle const& simfm = simParticles->at(trackId);  
      double trkrestm = simfm.endMomentum().e();
      if(i==0)_hEnergyat0->Fill(sqrt(mom.mag2() + trkrestm*trkrestm));
      if(i==1)_hEnergyat1->Fill(sqrt(mom.mag2() + trkrestm*trkrestm));
      
    } // end loop over hits.

    for ( size_t i=0; i<simParticles->size(); ++ i){
      SimParticle const& sim = simParticles->at(i);
      if(!sim.madeInG4())_hEnergysim->Fill(sim.startMomentum().e());
    }
  } // end FillHistograms




}  // end namespace mu2e

using mu2e::Histforpabs;
DEFINE_FWK_MODULE(Histforpabs);
