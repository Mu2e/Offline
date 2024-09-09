// Adapted from ReadVirtualDetector_module.cc
// Simplifies plotting for STM simulation studies. For a chosen virtual detector can generate a TTree with the hit time, pdgID, x and y position, and energy. Can summarise the data verbosely.
// Original author: Ivan Logashenko
// Adapted by: Pawel Plesniak

// stdlib includes
#include <cmath>
#include <iostream>
#include <string>

// exception handling
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"

// fhicl includes
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Utilities/InputTag.h"

// Offline includes
#include "Offline/ConditionsService/inc/ConditionsHandle.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/VirtualDetector.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/G4BeamlineInfo.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/Mu2eUtilities/inc/fromStrings.hh"

// CLHEP includes
#include "CLHEP/Units/SystemOfUnits.h"

// ROOT includes
#include "TH1F.h"
#include "TTree.h"
#include "TNtuple.h"

using namespace std;
namespace mu2e {
  class MakeTree : public art::EDAnalyzer {
    typedef SimParticleCollection::key_type key_type;
    art::ProductToken<StepPointMCCollection> StepPointMCsToken;
    art::ProductToken<SimParticleCollection> SimParticlemvToken;

    TTree* _ttree;
    float time = 0.0;
    int pdgId = 0, particle_count = 0;
    float x = 0.0, y = 0.0, z = 0.0;
    float mass = 0.0, E = 0.0;
    key_type trackId;

    bool GenerateTTree = true, GenerateDataSummary = true;
    StepPointMC::VolumeId_type VirtualDetectorID = 0;

    std::vector<int> pdgIdType;
    std::vector<int> pdgIdCount;
    std::vector<int>::iterator it;

  public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    struct Config {
      fhicl::Atom<art::InputTag> StepPointMCsTag{ Name("StepPointMCsTag"), Comment("Tag identifying the StepPointMCs")};
      fhicl::Atom<art::InputTag> SimParticlemvTag{ Name("SimParticlemvTag"), Comment("Tag identifying the SimParticlemv")};
      fhicl::Atom<bool> GenerateTTree{ Name("GenerateTTree"), Comment("Controls TTree generation")};
      fhicl::Atom<bool> GenerateDataSummary{ Name("GenerateDataSummary"), Comment("Controls dataset summary definition")};
      fhicl::OptionalAtom<StepPointMC::VolumeId_type> VirtualDetectorID{ Name("VirtualDetectorID"), Comment("Which virtual detector you check")};
    };
    using Parameters = art::EDAnalyzer::Table<Config>;

    explicit MakeTree(const Parameters& conf);
    void analyze(const art::Event& e);
    void endJob();

  };

  MakeTree::MakeTree(const Parameters& conf) :
    art::EDAnalyzer(conf),
    StepPointMCsToken(consumes<StepPointMCCollection>(conf().StepPointMCsTag())),
    SimParticlemvToken(consumes<SimParticleCollection>(conf().SimParticlemvTag())),
    GenerateTTree(conf().GenerateTTree()),
    GenerateDataSummary(conf().GenerateDataSummary())
  {
    if ( !( GenerateTTree || GenerateDataSummary ) ){
      throw cet::exception("CONFIGURATION") << "Either GenerateTTree or GenerateDataSummary need to be true.";
    }
    if (GenerateTTree){
      art::ServiceHandle<art::TFileService> tfs;
      _ttree = tfs->make<TTree>( "ttree", "Virtual Detectors ntuple");
      _ttree->Branch("time", &time, "time/F");
      _ttree->Branch("pdgId", &pdgId, "pdgId/I");
      _ttree->Branch("x", &x, "x/F");
      _ttree->Branch("y", &y, "y/F");
      _ttree->Branch("z", &z, "z/F");
      _ttree->Branch("E", &E, "E/F");
    }
    auto _VirtualDetectorID = conf().VirtualDetectorID();
    if(_VirtualDetectorID)VirtualDetectorID = *_VirtualDetectorID;
  };

  void MakeTree::analyze(const art::Event& event) {
    // Ask the event to give us a "handle" to the requested hits.
    auto const& StepPoints = event.getProduct(StepPointMCsToken);
    auto const& SimParticles = event.getProduct(SimParticlemvToken);
    GlobalConstantsHandle<ParticleDataList> pdt;

    // Loop over all VD hits.
    for (const StepPointMC& step : StepPoints){
      if((VirtualDetectorID != 0) && (step.volumeId() != VirtualDetectorID))
        continue;

      // Get the associated particle
      const SimParticle& particle = SimParticles.at(step.trackId());
      particle_count++;

      // Extract the parameters
      time = step.time();
      pdgId = particle.pdgId();
      x = step.position().x();
      y = step.position().y();
      z = step.position().z();
      mass = pdt->particle(particle.pdgId()).mass();
      E = std::sqrt(step.momentum().mag2()+mass*mass)-mass;

      // Generate the data summary
      if (GenerateDataSummary){
        it = std::find(pdgIdType.begin(), pdgIdType.end(), pdgId);
        if (it != pdgIdType.end())
          pdgIdCount[it - pdgIdType.begin()] = pdgIdCount[it - pdgIdType.begin()] + 1;
        else
          {
            pdgIdType.push_back(pdgId);
            pdgIdCount.push_back(1);
          };
      };
    };
    // Store all the data
    if (GenerateTTree)
      _ttree->Fill();
    return;
  };

  void MakeTree::endJob(){
    if (GenerateDataSummary){
      int c = pdgIdType.size();
      std::cout << "\n\nGeneral summary" << std::endl;
      std::cout << "Number of particles counted: " << particle_count << std::endl;
      std::cout << "Number of particle types: " << c << std::endl;
      for (int _i = 0; _i < c; _i++)
        std::cout << "PDGID " << pdgIdType[_i] << ": " << pdgIdCount[_i] << std::endl;
      std::cout << "\n\n" << std::endl;
    }
  };

}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::MakeTree)
