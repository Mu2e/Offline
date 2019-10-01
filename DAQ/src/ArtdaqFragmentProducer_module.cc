//
// EDAnalyzer module for reading in DataBlock data products from multiple
// root files and producing a binary output file compatible with the DTC
//
//
// Original author Tomo Miyashita
//
//

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>

#include <cassert>

// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Principal/Provenance.h"
#include "cetlib_except/exception.h"

// CLHEP includes.
#include "CLHEP/Random/RandGaussQ.h"

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "DAQDataProducts/inc/DataBlockCollection.hh"

#include "artdaq-core/Data/Fragment.hh"

#include "dtcInterfaceLib/DTC_Types.h"
#include "dtcInterfaceLib/DTC_Packets.h"
#include "dtcInterfaceLib/DTC.h"

#include "mu2e-artdaq-core/Overlays/mu2eFragment.hh"

#include "SeedService/inc/SeedService.hh"

#include <fstream>
#include <stdexcept>

namespace art {
  class ArtdaqFragmentProducer;
}

using art::ArtdaqFragmentProducer;

//--------------------------------------------------------------------

class art::ArtdaqFragmentProducer
 : public EDProducer
{

public:

  // --- C'tor/d'tor:
  explicit ArtdaqFragmentProducer(fhicl::ParameterSet const& pset);
  virtual ~ArtdaqFragmentProducer() { }

  // --- Production:
  virtual void produce( Event & );

  virtual void beginJob();
  virtual void endJob();

private:

  void flushBuffer();

  int _parseTRK;
  int _parseCAL;
  int _parseCRV;

  std::string _trkFragmentsTag;
  std::string _caloFragmentsTag;
  std::string _crvFragmentsTag;

  std::string _binaryDataBlockMakerModule;

  // Diagnostics level.
  int _diagLevel;

  // Limit on number of events for which there will be full printout.
  int _maxFullPrint;

}; // ArtdaqFragmentProducer

// ======================================================================

ArtdaqFragmentProducer::ArtdaqFragmentProducer(fhicl::ParameterSet const& pset):
  EDProducer(pset),
  _parseTRK(pset.get<int>("includeTracker",1)),
  _parseCAL(pset.get<int>("includeCalorimeter",1)),
  _parseCRV(pset.get<int>("includeCosmicRayVeto",0)),
//  _trkFragmentsTag{"daq:trk"},
//  _caloFragmentsTag{"daq:calo"},
//  _crvFragmentsTag{"daq:crv"},
  _trkFragmentsTag{"offlinetrk"},
  _caloFragmentsTag{"offlinecalo"},
  _crvFragmentsTag{"offlinecrv"},
  _binaryDataBlockMakerModule{"binaryOutput"},
  _diagLevel(pset.get<int>("diagLevel",0))
{
  produces<artdaq::Fragments>(_trkFragmentsTag);
  produces<artdaq::Fragments>(_caloFragmentsTag);
  produces<artdaq::Fragments>(_crvFragmentsTag);
}

void ArtdaqFragmentProducer::beginJob(){
  if ( _diagLevel > 0 ) {
    std::cout << "ArtdaqFragmentProducer Diaglevel: "
         << _diagLevel << " "
         << _maxFullPrint
	      << std::endl;
  }
}

void ArtdaqFragmentProducer::endJob(){
}

void ArtdaqFragmentProducer::produce(Event & evt) {

  if(_diagLevel>1) {
    std::cout << "ArtdaqFragmentProducer entering produce()" << std::endl;
  }

  bool gblDresult;
  art::Handle< std::vector<mu2e::DataBlock::adc_t> > binaryDataHandle;
  gblDresult = evt.getByLabel(_binaryDataBlockMakerModule, binaryDataHandle);

  if (!gblDresult) throw cet::exception("DATA") << " Missing binary DataBlock DataProduct";
  
  uint8_t const* current_ = reinterpret_cast<const uint8_t*>(&binaryDataHandle->at(0));

  auto resultTRK = std::make_unique<artdaq::Fragments>();
  auto resultCAL = std::make_unique<artdaq::Fragments>();
  auto resultCRV = std::make_unique<artdaq::Fragments>();
  
  // Each fragment has N super blocks--these super blocks are what
  // will be broken up into art::Events.  For a given super block,
  // there are M data blocks.
  
  // Increment through the data blocks of the current super block.
  auto const begin = current_;
  auto data = reinterpret_cast<char const*>(begin);

  char const* end = data;
  if(binaryDataHandle->size()>0) {
    end = reinterpret_cast<char const*>(current_ + sizeof(mu2e::DataBlock::adc_t)*binaryDataHandle->size());
  }

  size_t counter = 0;

  while (data < end) {

    if(_diagLevel>1) {
      std::cout << "ArtdaqFragmentProducer \tBeginning loop: " << counter << std::endl;
      counter++;
    }

    // Construct DTC_DataHeaderPacket to determine byte count of
    // current data block.
    DTCLib::DTC_DataPacket const dataPacket{data};
    DTCLib::DTC_DataHeaderPacket const headerPacket{dataPacket};

    auto const byteCount = headerPacket.GetByteCount();
    
    // Use byte count to calculate how many words the current data
    // block should occupy in the new fragment.
    auto const wordCount = byteCount / sizeof(artdaq::RawDataType);
    auto const packetSize = (byteCount % 8 == 0) ? wordCount : wordCount + 1;

    if (_parseTRK && headerPacket.GetSubsystem() == DTCLib::DTC_Subsystem::DTC_Subsystem_Tracker) {
      resultTRK->push_back(*artdaq::Fragment::dataFrag(headerPacket.GetTimestamp().GetTimestamp(true),
						       headerPacket.GetEVBMode(),  // Returns evbMode (see mu2e-docdb 4914)
						       reinterpret_cast<artdaq::RawDataType const*>(data), packetSize,
						       headerPacket.GetTimestamp().GetTimestamp(true))
			   .get());
    } else if(_parseCAL && headerPacket.GetSubsystem() == DTCLib::DTC_Subsystem::DTC_Subsystem_Calorimeter) {
      resultCAL->push_back(*artdaq::Fragment::dataFrag(headerPacket.GetTimestamp().GetTimestamp(true),
						       headerPacket.GetEVBMode(),  // Returns evbMode (see mu2e-docdb 4914)
						       reinterpret_cast<artdaq::RawDataType const*>(data), packetSize,
						       headerPacket.GetTimestamp().GetTimestamp(true))
			   .get());
    } else if(_parseCRV && headerPacket.GetSubsystem() == DTCLib::DTC_Subsystem::DTC_Subsystem_CRV) {
      resultCRV->push_back(*artdaq::Fragment::dataFrag(headerPacket.GetTimestamp().GetTimestamp(true),
						       headerPacket.GetEVBMode(),  // Returns evbMode (see mu2e-docdb 4914)
						       reinterpret_cast<artdaq::RawDataType const*>(data), packetSize,
						       headerPacket.GetTimestamp().GetTimestamp(true))
			   .get());
    }
    data += byteCount;
  }

  if(data != end) {
    throw art::Exception{art::errors::DataCorruption, "ArtdaqFragmentProducer::produce "}
                                                   << "The data pointer has shot past the 'end' pointer.";
  }
	
  // Store the artdaq::Fragment collections in the event
  evt.put(std::move(resultTRK),_trkFragmentsTag);
  evt.put(std::move(resultCAL),_caloFragmentsTag);
  evt.put(std::move(resultCRV),_crvFragmentsTag);

} // end of ::produce


// ======================================================================

DEFINE_ART_MODULE(ArtdaqFragmentProducer);

// ======================================================================
