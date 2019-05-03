//
// Find a product, modify it and rewrite it.
//
// Original author Rob Kutschke
//

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"

#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StepPointMC.hh"

//#include <iostream>
using namespace std;

namespace mu2e {

  class ModifyTrackSPM : public art::EDProducer {
  public:

    explicit ModifyTrackSPM(fhicl::ParameterSet const& pset);

    void produce( art::Event& e) override;

  private:
    art::InputTag tag_;

  };

  ModifyTrackSPM::ModifyTrackSPM(fhicl::ParameterSet const& pset ):
    art::EDProducer{pset},
    tag_(pset.get<std::string>("productTag")){
    //    std::cout << "DNBug:  have tag_:  " << tag_ << std::endl;
    //    std::cout << "DNBug:  it has instance:  " << tag_.instance() << std::endl;
    produces<std::vector<mu2e::StepPointMC> >( tag_.instance() );
    //    std::cout << "DNBug:  have said I will produce this tag." << std::endl;
  }

  void ModifyTrackSPM::produce(art::Event& event) {

    auto testp = event.getValidHandle<std::vector<mu2e::StepPointMC> >(tag_);

    //    std::cout << "DNBug have testp for event  " << std::endl;

    unique_ptr<std::vector<mu2e::StepPointMC> > prod = 
      std::make_unique<std::vector<mu2e::StepPointMC> >();

    //    std::cout << "DNBug have created the pointer prod  "  << std::endl;

    CLHEP::Hep3Vector offset(0.0,0.0,29.0);  // Fixing a problem in z0
    //    int count = 0;
    for ( auto oldTSPM : *testp ) {
      //      std::cout << "DNBug:  Step " << ++count << ":  ";
      mu2e::StepPointMC myTSPM(oldTSPM.simParticle(),
			       oldTSPM.volumeId(),
			       oldTSPM.totalEDep(),
			       oldTSPM.nonIonizingEDep(),
			       oldTSPM.time(),
			       oldTSPM.properTime(),
			       oldTSPM.position() + offset,
			       oldTSPM.momentum(),
			       oldTSPM.stepLength(),
			       oldTSPM.endProcessCode() );
      //      std::cout << ".  Made myTSPM." << std::endl;
      prod->push_back(myTSPM);
      //      std::cout << "DNBug:  Pushed it to the prod vector." << std::endl;
    }
    event.put(std::move(prod),tag_.instance());
    //    std::cout << "DNBug:  put it in the event" << std::endl;
  } // end ModifyTrackSPM::analyze

}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::ModifyTrackSPM)
