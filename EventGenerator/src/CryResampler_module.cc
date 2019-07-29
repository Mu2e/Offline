#include <iostream>
#include <string>
#include <cmath>
#include <memory>
#include <algorithm>

#include "cetlib_except/exception.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "SeedService/inc/SeedService.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "GeneralUtilities/inc/RSNTIO.hh"

#include "TTree.h"
#include "TFile.h"
#include "art_root_io/TFileService.h"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"

namespace mu2e {
  class CryResampler : public art::EDProducer {
    const mu2e::ParticleDataTable *pdt_;
    int verbosityLevel_;
    std::string inputFile_;
    std::string treeName_;
    std::string branchName_;

  struct recordedParticle {
    float x;
    float y;
    float z;
    float time;

    float px;
    float py;
    float pz;
    float pmag;
    float ek;

    float charge;
    int   pdgId;
    unsigned particleId;
    unsigned volumeCopyNumber;
  };

    recordedParticle particle_;
    TTree * particles_;

    public:
    explicit CryResampler(const fhicl::ParameterSet& pset);
    virtual void produce(art::Event& event);
  };

  //================================================================
  CryResampler::CryResampler(const fhicl::ParameterSet& pset) :
    art::EDProducer{pset}
    {
    verbosityLevel_ = pset.get<int>("verbosityLevel", 0);
    pdt_ = (&*GlobalConstantsHandle<ParticleDataTable>());
    inputFile_ = pset.get<std::string>("inputFile");
    treeName_ = pset.get<std::string>("treeName", "crvDump/nt");
    branchName_ = pset.get<std::string>("branchName", "hits");

    art::ServiceHandle<art::TFileService> tfs;
    const std::string resolvedFileName = ConfigFileLookupPolicy()(inputFile_);
    TFile *input = tfs->make<TFile>(resolvedFileName.c_str(), "READ");
    particles_ = (TTree*)input->Get(treeName_.c_str());

    if(!particles_) {
      throw cet::exception("BADINPUT")<<"CryResampler: Could not get tree \""
        << treeName_ <<"\" from file \""<<input->GetName() <<"\"\n";
    }

    particles_->SetBranchAddress(branchName_.c_str(), &particle_);


    std::cout << "Using " << particles_->GetEntries() << " records from " <<
      inputFile_ << std::endl;


    produces<GenParticleCollection>();
  }

  //================================================================
  void CryResampler::produce(art::Event& event) {
    std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);

    particles_->GetEntry(event.id().event());

    const CLHEP::Hep3Vector p3(particle_.px, particle_.py, particle_.pz);
    const double mass = pdt_->particle(particle_.pdgId).ref().mass().value();
    const double energy = std::sqrt(std::pow(mass,2) + p3.mag2());
    CLHEP::HepLorentzVector fourmom(p3, energy);
    const CLHEP::Hep3Vector pos(particle_.x, particle_.y, particle_.z);
    // output->emplace_back(PDGCode::type(particle_.pdgId),
        // GenId::cosmicCRY,
        // pos,
        // fourmom,
        // particle_.time);
    output->push_back(GenParticle(static_cast<PDGCode::type>(particle_.pdgId),
              GenId::cosmicCRY, pos, fourmom, 0.));

    if (verbosityLevel_ > 1) {
      std::cout << event.id().event() << ": "
        << particle_.pdgId << ", "
        << p3 << ", " << mass << ", " << energy
        << std::endl;
    }

    event.put(std::move(output));
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::CryResampler);
