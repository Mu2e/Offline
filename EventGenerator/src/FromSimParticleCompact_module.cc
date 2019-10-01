// Continue with compact sim particles (a Root Tree) from previous jobs
//
// Zhengyun You, 2013-12-01

#include <iostream>
#include <string>
#include <cmath>
#include <memory>
#include <algorithm>

#include "cetlib_except/exception.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"

#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "SeedService/inc/SeedService.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Mu2eUtilities/inc/CzarneckiSpectrum.hh"
#include "Mu2eUtilities/inc/SimpleSpectrum.hh"
#include "Mu2eUtilities/inc/EjectedProtonSpectrum.hh"
#include "Mu2eUtilities/inc/BinnedSpectrum.hh"
#include "Mu2eUtilities/inc/Table.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "ProductionTargetGeom/inc/ProductionTarget.hh"

#include "TTree.h"
#include "TFile.h"

namespace mu2e {

  //================================================================
  namespace {
    struct VDHit {
      float x;
      float y;
      float z;
      float time;
      float px;
      float py;
      float pz;
      int   pdgId;

      VDHit() : x(std::numeric_limits<double>::quiet_NaN())
              , y(std::numeric_limits<double>::quiet_NaN())
              , z(std::numeric_limits<double>::quiet_NaN())
              , time(std::numeric_limits<double>::quiet_NaN())
              , px(std::numeric_limits<double>::quiet_NaN())
              , py(std::numeric_limits<double>::quiet_NaN())
              , pz(std::numeric_limits<double>::quiet_NaN())
              , pdgId(0)
      {}

    }; // struct VDHit
  } // namespace {}

  //================================================================
  class FromSimParticleCompact : public art::EDProducer {
    std::vector<std::string> inputFiles_;
    std::string treeName_;

    PDGCode::type pdgId_;
    int verbosityLevel_;

    std::vector<VDHit> vdhits_;
    void loadInputFiles();
    void addVDHits(TFile* infile);

    bool rotateTarget_;

    double phiMin_;   // rotate an angle between phiMin_ and phiMax_
    double phiMax_;

    art::RandomNumberGenerator::base_engine_t& eng_;
    CLHEP::RandFlat randFlat_;

    CLHEP::HepRotation targetRotation_; // rotates target frame to Mu2e frame
    CLHEP::Hep3Vector  targetCenter_;

    bool firstEvent_;

  public:
    explicit FromSimParticleCompact(const fhicl::ParameterSet& pset);
    virtual void produce(art::Event& event);
    virtual void beginRun(art::Run& run);
  };

  //================================================================
  FromSimParticleCompact::FromSimParticleCompact(const fhicl::ParameterSet& pset)
    : EDProducer{pset}
    , inputFiles_(pset.get<std::vector<std::string> >("inputFiles"))
    , treeName_(pset.get<std::string>("treeName"))
    , verbosityLevel_(pset.get<int>("verbosityLevel", 0))
    , rotateTarget_(pset.get<bool>("rotateTarget", false))
    , phiMin_(pset.get<double>("phiMin", 0.))
    , phiMax_(pset.get<double>("phiMax", 360.))
    , eng_(createEngine(art::ServiceHandle<SeedService>()->getSeed()))
    , randFlat_(eng_)
    , firstEvent_(true)
  {
    produces<mu2e::GenParticleCollection>();
    if(inputFiles_.empty()) {
      throw cet::exception("BADCONFIG")<<"Error: no inputFiles";
    }

    if(verbosityLevel_ > 0) {
      std::cout<<"FromSimParticleCompact: producing particle "
               <<"..."
               <<std::endl;
    }

    phiMin_   *= CLHEP::deg;
    phiMax_   *= CLHEP::deg;
  }

  //================================================================
  void FromSimParticleCompact::beginRun(art::Run& run) {
    loadInputFiles();

    if(verbosityLevel_ > 0) {
      std::cout<<"FromSimParticleCompact inputs: num input particles = "
               <<vdhits_.size()
               <<std::endl;
    }
  }

  //================================================================
  void FromSimParticleCompact::loadInputFiles() {
    art::ServiceHandle<art::TFileService> tfs;

    for(const auto& fn : inputFiles_) {
      const std::string resolvedFileName = ConfigFileLookupPolicy()(fn);

      if(verbosityLevel_ > 0) {
        std::cout<<"FromSimParticleCompact: reading input file "<<resolvedFileName<<std::endl;
      }

      TFile *infile = tfs->make<TFile>(resolvedFileName.c_str(), "READ");
      addVDHits(infile);
    }

  } // loadInputFiles()

  //================================================================
  void FromSimParticleCompact::addVDHits(TFile* infile) {
    TTree *nt = dynamic_cast<TTree*>(infile->Get(treeName_.c_str()));
    if(!nt) {
      throw cet::exception("BADINPUT")<<"Could not get tree \""<<treeName_
                                      <<"\" from file \""<<infile->GetName()
                                      <<"\"\n";
    }

    const Long64_t nTreeEntries = nt->GetEntries();
    if(verbosityLevel_>0) {
      std::cout<<"FromSimParticleCompact: nTreeEntries for particles = "<<nTreeEntries<<std::endl;
    }

    VDHit hit;
    TBranch *bhit = nt->GetBranch("particles");
    bhit->SetAddress(&hit);

    vdhits_.reserve(vdhits_.size() + nTreeEntries);
    for(Long64_t i=0; i<nTreeEntries; ++i) {
      bhit->GetEntry(i);
      vdhits_.emplace_back(hit);
    }
  }

  //================================================================
  void FromSimParticleCompact::produce(art::Event& event) {

    if (firstEvent_) {
      GeomHandle<ProductionTarget> target;
      targetRotation_ = target->protonBeamRotation();
      targetCenter_ = target->position();
      if(verbosityLevel_ > 1) {
        std::cout << "targetCenter   " << targetCenter_ << std::endl;
        std::cout << "targetRotation " << targetRotation_ << std::endl;
      }
      firstEvent_ = false;
    }

    std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);

    const VDHit& hit = vdhits_.at(randFlat_.fireInt(vdhits_.size()));

    const CLHEP::Hep3Vector pos(hit.x, hit.y, hit.z);
    const CLHEP::Hep3Vector mom(hit.px, hit.py, hit.pz);

    CLHEP::Hep3Vector mom_new = mom;
    CLHEP::Hep3Vector pos_new = pos;

    if (rotateTarget_) {
      if(verbosityLevel_ > 1) std::cout << "pos " << pos << " mom " << mom << std::endl;

      CLHEP::Hep3Vector mom_tgt = targetRotation_.inverse()*mom;
      CLHEP::Hep3Vector pos_tgt = targetRotation_.inverse()*(pos-targetCenter_);
      if(verbosityLevel_ > 1) std::cout << "pos_tgt " << pos_tgt << " mom_tgt " << mom_tgt << std::endl;

      // rotate by an random angle around phi and get new p direction in target frame
      double rotatePhi = randFlat_.fire()*(phiMax_-phiMin_);
      if(verbosityLevel_ > 1) std::cout << "rotatePhi " << rotatePhi << std::endl;

      mom_tgt.rotateZ(rotatePhi);
      pos_tgt.rotateZ(rotatePhi);
      if(verbosityLevel_ > 1) std::cout << "pos_tgt_rot " << pos_tgt << " mom_tgt_rot " << mom_tgt << std::endl;

      // rotate around y back to mu2e frame
      mom_new = targetRotation_*mom_tgt;
      pos_new = targetRotation_*pos_tgt+targetCenter_;
      if(verbosityLevel_ > 1) std::cout << "pos_rot " << pos_new << " mom_rot " << mom_new << std::endl;
    }

    PDGCode::type pdgId = static_cast<PDGCode::type>(hit.pdgId);
    const double mass = GlobalConstantsHandle<ParticleDataTable>()->particle(pdgId).ref().mass().value();
    const double energy = sqrt(mom_new.mag2() + mass*mass);
    CLHEP::HepLorentzVector fourmom(mom_new, energy);
    const double time = hit.time;

    GenParticle outGen(pdgId, GenId::fromSimParticleCompact, pos_new, fourmom, time);
    output->push_back(outGen);

    event.put(std::move(output));
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::FromSimParticleCompact);
