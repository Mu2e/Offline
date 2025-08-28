// Generate an example lumi stream information for DAQ development

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "Offline/RecoDataProducts/inc/IntensityInfoCalo.hh"
#include "Offline/RecoDataProducts/inc/IntensityInfoTimeCluster.hh"
#include "Offline/RecoDataProducts/inc/IntensityInfoTrackerHits.hh"

#include <iostream>

namespace mu2e
{

  class DummyLumiInfoProducer : public art::EDProducer
  {
  public:
    struct Config
    {
      fhicl::Atom<int> diagLevel{fhicl::Name("diagLevel"), fhicl::Comment("diagnostic Level"), 0};
      fhicl::Atom<int> simMode{fhicl::Name("simMode"), fhicl::Comment("Simulation mode: 0 = all zeros; 1 = non-zero values"), 0};
    };

    explicit DummyLumiInfoProducer(const art::EDProducer::Table<Config>& config);
    virtual void produce(art::Event& event) override;

  private:
    int   _diagLevel;
    int   _simMode;
  };

  DummyLumiInfoProducer::DummyLumiInfoProducer(const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config}
    , _diagLevel(config().diagLevel())
    , _simMode(config().simMode())
  {
    produces<mu2e::IntensityInfoCalo       >();
    produces<mu2e::IntensityInfoTimeCluster>();
    produces<mu2e::IntensityInfoTrackerHits>();
  }

  void DummyLumiInfoProducer::produce(art::Event& event)
  {
    const art::EventNumber_t  eventNumber  = event.event ();
    const art::SubRunNumber_t subrunNumber = event.subRun();
    const art::RunNumber_t    runNumber    = event.run   ();

    if(_diagLevel > 1) std::cout << "DummyLumiInfoProducer::" << __func__ << ": Event " << runNumber << ":" << subrunNumber << ":" << eventNumber << std::endl;

    //---------------------------------------
    // Create example data

    auto caloInfo        = std::unique_ptr<mu2e::IntensityInfoCalo       >(new mu2e::IntensityInfoCalo       );
    auto timeClusterInfo = std::unique_ptr<mu2e::IntensityInfoTimeCluster>(new mu2e::IntensityInfoTimeCluster);
    auto trackerInfo     = std::unique_ptr<mu2e::IntensityInfoTrackerHits>(new mu2e::IntensityInfoTrackerHits);

    // Assign non-zero entries for harder compression
    if(_simMode == 1) {
      constexpr int prime_1(17), prime_2(251), prime_3(503), prime_4(1523); //for distributing values somewhat evenly
      caloInfo->setNCaloHitsD0((eventNumber * prime_2) % (prime_4));
      caloInfo->setNCaloHitsD1((eventNumber * prime_2) % (prime_3));
      caloInfo->setCaloEnergy ((eventNumber * prime_1) % (prime_2));
      const int nCaphriHits = (eventNumber * prime_3) % (prime_1);
      constexpr int max_caphri_index = 4; // only four CAPHRI crystals
      std::vector<unsigned short> caphriHits;
      for(int ihit = 0; ihit < nCaphriHits; ++ihit) {
        const int index = CaloConst::_caphriId[(eventNumber * prime_1) % max_caphri_index];
        const double energy = (eventNumber * prime_2) % prime_1;
        caloInfo->addCaphriHit(energy, index);
      }

      timeClusterInfo->setNProtonTCs((eventNumber * prime_1) % (prime_2));

      trackerInfo->setNTrackerHits((eventNumber * prime_1) % (prime_4));
    }

    // Add the intensity info to the event
    event.put(std::move(caloInfo));
    event.put(std::move(timeClusterInfo));
    event.put(std::move(trackerInfo));

  } //produce

} //namespace mu2e

using mu2e::DummyLumiInfoProducer;
DEFINE_ART_MODULE(DummyLumiInfoProducer)
