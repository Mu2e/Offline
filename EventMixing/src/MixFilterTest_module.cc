#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Modules/MixFilter.h"
#include "art/Framework/IO/ProductMix/MixHelper.h"
#include "art/Framework/Core/PtrRemapper.h"
#include "art/Persistency/Common/CollectionUtilities.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "cetlib/map_vector.h"
#include "canvas/Utilities/InputTag.h"

#include "boost/noncopyable.hpp"

namespace arttest {
  class MixFilterTestDetail;
  typedef art::MixFilter<MixFilterTestDetail> MixFilterTest;
}

class arttest::MixFilterTestDetail : private boost::noncopyable {
public:
  // Constructor is responsible for registering mix operations with
  // MixHelper::declareMixOp() and bookkeeping products with
  // MixHelperproduces().
  MixFilterTestDetail(fhicl::ParameterSet const &p,
                      art::MixHelper &helper);    

  // Return the number of secondaries to read this time. Declare const
  // if you don't plan to change your class' state.
  size_t nSecondaries() const;

  // Optional processEventIDs(): after the generation of the event
  // sequence, this function will be called if it exists to provide the
  // sequence of EventIDs.
  void processEventIDs(art::EventIDSequence const &seq);

  // Optional.finalizeEvent(): (eg) put bookkeping products in event. Do
  // *not* place mix products into the event: this will already have
  // been done for you.
  void finalizeEvent(art::Event &t);


  bool
  aggregateDoubleCollection(std::vector<std::vector<double> const *> const &in,
                      std::vector<double> &out,
                      art::PtrRemapper const &);

private:
  size_t nSecondaries_;
  bool testRemapper_;
  std::vector<size_t> doubleVectorOffsets_;
  std::unique_ptr<art::EventIDSequence> eIDs_;
  bool startEvent_called_;
  bool processEventIDs_called_;
};

arttest::MixFilterTestDetail::
MixFilterTestDetail(fhicl::ParameterSet const &p,
                      art::MixHelper &helper)
  :
  nSecondaries_(p.get<size_t>("numSecondaries", 1)),
  testRemapper_(p.get<bool>("testRemapper", 1)),
  doubleVectorOffsets_(),
  eIDs_(),
  startEvent_called_(false),
  processEventIDs_called_(false)
{
  std::string mixProducerLabel(p.get<std::string>("mixProducerLabel",
                                                  "mixProducer"));

  helper.produces<std::string>(); // "Bookkeeping"
  helper.produces<art::EventIDSequence>(); // "Bookkeeping"

  helper.declareMixOp
    (art::InputTag(mixProducerLabel, "doubleCollectionLabel"),
     &MixFilterTestDetail::aggregateDoubleCollection, *this);

}


size_t
arttest::MixFilterTestDetail::
nSecondaries() const {
  return nSecondaries_;
}

void
arttest::MixFilterTestDetail::
processEventIDs(art::EventIDSequence const &seq) {
  processEventIDs_called_ = true;
  eIDs_.reset(new art::EventIDSequence(seq));
}

void
arttest::MixFilterTestDetail::
finalizeEvent(art::Event &e) {
  e.put(std::unique_ptr<std::string>(new std::string("BlahBlahBlah")));
  e.put(std::move(eIDs_));

  assert(startEvent_called_);
  assert(processEventIDs_called_);
  startEvent_called_ = false;
  processEventIDs_called_ = false;
}

bool
arttest::MixFilterTestDetail::
aggregateDoubleCollection(std::vector<std::vector<double> const *> const &in,
                    std::vector<double> &out,
                    art::PtrRemapper const &) {
  art::flattenCollections(in, out, doubleVectorOffsets_);
  return true; //  Always want product in event.
}


DEFINE_ART_MODULE(arttest::MixFilterTest);

