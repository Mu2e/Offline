// Create G4 excited ions like those we have to save and restore
// sometimes in Mu2eG4 multistage jobs.
// This is for testing the software functionality.
//
// Andrei Gaponenko, 2020

#include <memory>

#include "cetlib_except/exception.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/TupleAs.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"

#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/PhysicalVolumeInfoMultiCollection.hh"

#include "Geant4/G4IonTable.hh"

namespace mu2e {

  //================================================================
  class IonProducer: public art::EDProducer {
    PDGCode::type pdgId_;
    double excitationEnergy_;
    short int floatLevel_;
    CLHEP::Hep3Vector pos_;
    CLHEP::Hep3Vector mom_;
  public:

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Atom<int> pdgId { Name("pdgId"), Comment("PDG ID of ion to create.") };
      fhicl::Atom<double> excitationEnergy { Name("excitationEnergy"), Comment("Ion excitation energy") };
      fhicl::Atom<short> floatLevel { Name("floatLevel"), Comment("float level base index")};

      fhicl::TupleAs<CLHEP::Hep3Vector(double, double, double)> position{ Name("position"), Comment("Where to place the ion") };
      fhicl::TupleAs<CLHEP::Hep3Vector(double, double, double)> momentum{ Name("momentum"), Comment("Three-momentum of the ion") };
    };

    using Parameters = art::EDProducer::Table<Config>;
    explicit IonProducer(const Parameters& conf);

    virtual void produce(art::Event& event) override;
    virtual void beginSubRun(art::SubRun& sr) override;
  };

  //================================================================
  IonProducer::IonProducer(const Parameters& conf)
    : EDProducer{conf}
    , pdgId_(PDGCode::type(conf().pdgId()))
    , excitationEnergy_(conf().excitationEnergy())
    , floatLevel_(conf().floatLevel())
    , pos_(conf().position())
    , mom_(conf().momentum())
  {
    produces<mu2e::SimParticleCollection>();
    produces<PhysicalVolumeInfoMultiCollection,art::InSubRun>();

    if (pdgId_ <= PDGCode::G4Threshold) {
      throw cet::exception("BADCONFIG") << "IonProducer: pdgId = "<<" is to small for an ion " << std::endl;
    }

  }

  //================================================================
  void IonProducer::beginSubRun(art::SubRun& sr) {
    using Collection_t = PhysicalVolumeInfoMultiCollection;
    auto mvi = std::make_unique<Collection_t>();

    // Create an entry for one stage.  We do not need volume info
    // content for the ion test use, but we still need a stage count.
    mvi->resize(1);

    sr.put(std::move(mvi), art::fullSubRun());
  }

  //================================================================
  void IonProducer::produce(art::Event& event) {

    auto particles = std::make_unique<SimParticleCollection>();

    // get A & A based on pdgId
    G4int ZZ,AA, lvl;
    G4double EE(0.0);
    if(!G4IonTable::GetNucleusByEncoding(pdgId_,ZZ,AA,EE,lvl)) {
      throw cet::exception("RUNTIME") << "IonProducer: G4IonTable did not resolve pdgId = "<<pdgId_<< std::endl;
    }

    const double mass = G4IonTable::GetIonTable()->GetIonMass(ZZ,AA);

    CLHEP::HepLorentzVector fourmom(mom_, mass);

    SimParticle::IonDetail ion(excitationEnergy_, floatLevel_);

    cet::map_vector_key one{1};
    SimParticle sp(one,  // index
                   0,    // simStage
                   art::Ptr<SimParticle>(), // parent
                   pdgId_,
                   art::Ptr<GenParticle>(),
                   pos_,
                   fourmom,
                   0, // startGlobalTime
                   0., // properTime
                   0, // start volume index
                   1, // G4 status
                   ProcessCode::mu2ePrimary,
                   ion
                   );

    sp.addEndInfo(pos_,
                  fourmom,
                  0., // time
                  0., //proper time
                  0, // volume index
                  1, // G4 status
                  ProcessCode::mu2ePrimary,
                  fourmom.e() - mass, // preLastStepKE
                  0, // nsteps
                  0. // track length
                  );

    particles->insert(std::make_pair(one, sp));

    event.put(std::move(particles));
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::IonProducer)
