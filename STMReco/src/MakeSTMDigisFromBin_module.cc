//
// Create STMDigis from STMSteps
//
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"

#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include <utility>
// root
#include "TH1F.h"
#include "TTree.h"

#include "Offline/MCDataProducts/inc/STMStep.hh"
#include "Offline/RecoDataProducts/inc/STMDigi.hh"

#include <bitset>

using namespace std;
using CLHEP::Hep3Vector;
namespace mu2e {

  struct Test {
    uint16_t trigNum;
    uint16_t trigType;
    uint32_t trigTime;
    uint16_t trigTimeOffset;
    uint16_t baselineMean;
    uint16_t baselineRMS;
    uint16_t nDrop;
    uint16_t ADC0;
    uint16_t ADC1;
    uint16_t ADC2;
    uint16_t ADC3;

    friend std::ostream& operator<<(ostream& os, const Test& test) {
      os << "Trigger #" << test.trigNum << std::endl;
      os << "\tTrigger Type: " << test.trigType << std::endl;
      os << "\tTrigger Time: " << test.trigTime << std::endl;
      os << "\tTrigger Time Offset: " << test.trigTimeOffset << std::endl;
      os << "\tBaseline Mean: " << test.baselineMean << std::endl;
      os << "\tBaseline RMS: " << test.baselineRMS << std::endl;
      os << "\tNo. of Dropped Packets: " << test.nDrop << std::endl;
      os << "\tADCs: " << test.ADC0 << " " << test.ADC1 << " " << test.ADC2 << " " << test.ADC3 << std::endl;
      return os;
    }
  };

  class MakeSTMDigisFromBin : public art::EDProducer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
	fhicl::Atom<std::string> binFileName{ Name("binFileName"), Comment("Binary file name")};
      };
      using Parameters = art::EDProducer::Table<Config>;
      explicit MakeSTMDigisFromBin(const Parameters& conf);

    private:
      void produce(art::Event& e) override;

    std::string _binFileName;
  };

  MakeSTMDigisFromBin::MakeSTMDigisFromBin(const Parameters& config )  : 
    art::EDProducer{config},
    _binFileName(config().binFileName())
  {
    produces<STMDigiCollection>();
  }

    void MakeSTMDigisFromBin::produce(art::Event& event) {
      // create output
      unique_ptr<STMDigiCollection> outputSTMDigis(new STMDigiCollection);

      std::ifstream binfile;
      binfile.open(_binFileName.c_str(), ios::binary);
      if (!binfile.is_open()) {
	throw cet::exception("MakeSTMDigisFromBin") << "A problem opening binary file " << _binFileName << std::endl;
      }
      Test test[1];
      while (binfile.read((char *) &test[0], sizeof(Test))) {
	std::cout << test[0] << std::endl;
	if (!binfile) {
	  throw cet::exception("MakeSTMDigisFromBin") << "A problem reading binary file " << _binFileName << std::endl;
	}
	STMDigi stm_digi(test[0].trigTimeOffset, test[0].ADC0);
	outputSTMDigis->push_back(stm_digi);
      }
      binfile.close();

      event.put(std::move(outputSTMDigis));
  }
}

DEFINE_ART_MODULE(mu2e::MakeSTMDigisFromBin)
