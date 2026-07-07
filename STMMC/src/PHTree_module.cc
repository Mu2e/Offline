// Adapted from ReadVirtualDetector_module.cc
// For STMPHDigis, generates a TTree with the time and energy
// Original author: Ivan Logashenko
// Adapted by: Pawel Plesniak

// stdlib includes
#include <cmath>
#include <iostream>

// art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"

// exception handling
#include "cetlib_except/exception.h"

// fhicl includes
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"

// message handling
#include "messagefacility/MessageLogger/MessageLogger.h"

// Offline includes
#include "Offline/RecoDataProducts/inc/STMPHDigi.hh"

// ROOT includes
#include "art_root_io/TFileService.h"
#include "TTree.h"


typedef cet::map_vector_key key_type;
typedef unsigned long VolumeId_type;

namespace mu2e {
  class PHTree : public art::EDAnalyzer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<art::InputTag> STMPHDigiTag{Name("STMPHDigiTag"), Comment("Tag identifying the PH Digis")};
        fhicl::OptionalAtom<double> EnergyCalib{ Name("EnergyCalib"), Comment("Controls whether to make the energy TTrees with units of energy or in ADC values. If 0, will leave as ADC value, otherwise will multiply by this calibration to generate the energy.")};
      };
      using Parameters = art::EDAnalyzer::Table<Config>;
      explicit PHTree(const Parameters& conf);
      void analyze(const art::Event& e);
      void endJob();
    private:
      art::ProductToken<STMPHDigiCollection> STMPHDigiToken;
      int eventCounter = 0, digiCounter = 0;
      TTree* ttree;
      uint32_t time = 0;
      double E = 0, EnergyCalib = 0;
  };

  PHTree::PHTree(const Parameters& conf) :
    art::EDAnalyzer(conf),
    STMPHDigiToken(consumes<STMPHDigiCollection>(conf().STMPHDigiTag())) {
      EnergyCalib = conf().EnergyCalib() ? *(conf().EnergyCalib()) : 1.0;
      art::ServiceHandle<art::TFileService> tfs;
      ttree = tfs->make<TTree>("ttree", "PH ttree");
      ttree->Branch("time", &time, "time/i"); // ns
      ttree->Branch("E", &E, "E/D"); // keV
    };

  void PHTree::analyze(const art::Event& event) {
    // Get the data products from the event
    auto const& PHDigis = event.getProduct(STMPHDigiToken);
    if (PHDigis.empty())
      return;
    eventCounter++;
    // Loop over all VD hits
    for (const STMPHDigi& digi : PHDigis) {
      // Extract the parameters
      time = digi.time();
      E = digi.energy() * EnergyCalib;
      if (E < 0) {
        std::cout << "Energy: " << E << std::endl;
        throw cet::exception("LogicError", "Energy is negative");
      };
      ttree->Fill();
      digiCounter++;
    };
    return;
  };

  void PHTree::endJob() {
    mf::LogInfo log("PH tree summary");
    log << "========= PH tree summary =========\n";
    log << std::left << std::setw(25) << "\tProcessed events: " << eventCounter << "\n";
    log << std::left << std::setw(25) << "\tProcessed digis: "  << digiCounter << "\n";
    log << "===================================\n";
  };
}; // end namespace mu2e

DEFINE_ART_MODULE(mu2e::PHTree)
