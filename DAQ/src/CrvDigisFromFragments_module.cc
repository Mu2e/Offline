// ======================================================================
//
// CrvDigisFromFragments_plugin:  Add CRV data products to the event
//
// ======================================================================

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Principal/Handle.h"
#include "mu2e-artdaq-core/Overlays/ArtFragmentReader.hh"

#include <artdaq-core/Data/Fragment.hh>
#include "RecoDataProducts/inc/StrawDigiCollection.hh"
#include "RecoDataProducts/inc/CaloDigiCollection.hh"
#include "RecoDataProducts/inc/CrvDigiCollection.hh"

#include <iostream>

#include <string>

#include <memory>

namespace art {
  class CrvDigisFromFragments;
}

using art::CrvDigisFromFragments;

// ======================================================================

class art::CrvDigisFromFragments
  : public EDProducer
{

public:

  using EventNumber_t = art::EventNumber_t;
  using adc_t = mu2e::ArtFragmentReader::adc_t;
  
  // --- C'tor/d'tor:
  explicit  CrvDigisFromFragments(fhicl::ParameterSet const& pset);
  virtual  ~CrvDigisFromFragments()  { }

  // --- Production:
  virtual void produce( Event & );

private:
  int   diagLevel_;

  int   parseCRV_;

  art::InputTag crvFragmentsTag_;

};  // CrvDigisFromFragments

// ======================================================================

CrvDigisFromFragments::CrvDigisFromFragments(fhicl::ParameterSet const& pset)
  : EDProducer{pset}
  , diagLevel_(pset.get<int>("diagLevel",0))
  , parseCRV_(pset.get<int>("parseCRV",1))
  , crvFragmentsTag_(pset.get<art::InputTag>("crvTag","daq:crv"))
{
  produces<EventNumber_t>(); 
  produces<mu2e::CrvDigiCollection>();
}

// ----------------------------------------------------------------------

void
  CrvDigisFromFragments::
  produce( Event & event )
{

  art::EventNumber_t eventNumber = event.event();

  auto crvFragments = event.getValidHandle<artdaq::Fragments>(crvFragmentsTag_);
  size_t numCrvFrags = crvFragments->size();

  if( diagLevel_ > 1 ) {
    std::cout << std::dec << "Producer: Run " << event.run() << ", subrun " << event.subRun()
	      << ", event " << eventNumber << " has " << std::endl;
    std::cout << crvFragments->size() << " CRV fragments." << std::endl;

    size_t totalSize = 0;
    for(size_t idx = 0; idx < crvFragments->size(); ++idx) {
      auto size = ((*crvFragments)[idx]).size() * sizeof(artdaq::RawDataType);
      totalSize += size;
      //      std::cout << "\tCRV Fragment " << idx << " has size " << size << std::endl;
    }
    
    std::cout << "\tTotal Size: " << (int)totalSize << " bytes." << std::endl;  
  }

  // Collection of CaloDigis for the event
  std::unique_ptr<mu2e::CrvDigiCollection> crv_digis(new mu2e::CrvDigiCollection);

  // Loop over the CRV fragments
  for (size_t idx = 0; idx < numCrvFrags; ++idx) {

    const auto& fragment((*crvFragments)[idx]);
    mu2e::ArtFragmentReader cc(fragment);
    
    if( diagLevel_ > 1 ) {
      std::cout << std::endl;
      std::cout << "ArtFragmentReader: ";
      std::cout << "\tBlock Count: " << std::dec << cc.block_count() << std::endl;
      std::cout << "\tByte Count: " << cc.byte_count() << std::endl;
      std::cout << std::endl;
      std::cout << "\t" << "====== Example Block Sizes ======" << std::endl;
      for(size_t i=0; i<10; i++) {
	if(i <cc.block_count()) {
	  std::cout << "\t" << i << "\t" << cc.blockIndexBytes(i) << "\t" << cc.blockSizeBytes(i) << std::endl;
	}
      }
      std::cout << "\t" << "=========================" << std::endl;
    }
    
    for(size_t curBlockIdx=0; curBlockIdx<cc.block_count(); curBlockIdx++) {

      size_t blockStartBytes = cc.blockIndexBytes(curBlockIdx);
      size_t blockEndBytes = cc.blockEndBytes(curBlockIdx);

      if( diagLevel_ > 1 ) {
	std::cout << "BLOCKSTARTEND: " << blockStartBytes << " " << blockEndBytes << " " << cc.blockSizeBytes(curBlockIdx)<< std::endl;
	std::cout << "IndexComparison: " << cc.blockIndexBytes(0)+16*(0+3*curBlockIdx) << "\t";
	std::cout                        << cc.blockIndexBytes(curBlockIdx)+16*(0+3*0) << std::endl;
      }

      adc_t const *pos = reinterpret_cast<adc_t const *>(cc.dataAtBytes(blockStartBytes));

      if( diagLevel_ > 1 ) {
	// Print binary contents the first 3 packets starting at the current position
	// In the case of the tracker simulation, this will be the whole tracker
	// DataBlock. In the case of the calorimeter, the number of data packets
	// following the header packet is variable.
	cc.printPacketAtByte(cc.blockIndexBytes(0)+16*(0+3*curBlockIdx));
	cc.printPacketAtByte(cc.blockIndexBytes(0)+16*(1+3*curBlockIdx));
	cc.printPacketAtByte(cc.blockIndexBytes(0)+16*(2+3*curBlockIdx));
	
	// Print out decimal values of 16 bit chunks of packet data
	for(int i=7; i>=0; i--) {
	  std::cout << (adc_t) *(pos+i);
	  std::cout << " ";
	}
	std::cout << std::endl;
      }	    

      adc_t rocID = cc.DBH_ROCID(pos);
      adc_t valid = cc.DBH_Valid(pos);
      adc_t packetCount = cc.DBH_PacketCount(pos);
	    
      uint32_t timestampLow    = cc.DBH_TimestampLow(pos);
      uint32_t timestampMedium = cc.DBH_TimestampMedium(pos);
      size_t timestamp = timestampLow | (timestampMedium<<16);
      
      adc_t EVBMode = cc.DBH_EVBMode(pos);
      adc_t sysID = cc.DBH_SubsystemID(pos);
      adc_t dtcID = cc.DBH_DTCID(pos);

      eventNumber = timestamp;
      
      if(sysID!=2) {
	throw cet::exception("DATA") << " CRV packet does not have system ID 2";
      }

      // Parse phyiscs information from the CRV packets
      if(packetCount>0 && parseCRV_>0) {

	for(size_t i=0; i<cc.DBVR_NumHits(pos); i++) {

	  // Fill the CrvDigiCollection
	  // CrvDigi(const std::array<unsigned int, NSamples> &ADCs, unsigned int startTDC,
	  //         mu2e::CRSScintillatorBarIndex scintillatorBarIndex, int SiPMNumber) :
	  crv_digis->emplace_back(cc.DBV_ADCs(pos,i),
				  cc.DBV_startTDC(pos,i),
				  mu2e::CRSScintillatorBarIndex(cc.DBV_sipmID(pos,i)/4),
				  cc.DBV_sipmID(pos,i)%4);
	}

	if( diagLevel_ > 1 ) {

	  for(size_t i=0; i<cc.DBVR_NumHits(pos); i++) {
	    std::cout << "MAKEDIGI: " << cc.DBV_sipmID(pos,i)%4 << " " << cc.DBV_sipmID(pos,i)/4 << " " << cc.DBV_startTDC(pos,i)
		      << " " << cc.DBVR_NumHits(pos) << " ";
	    
	    for(size_t j=0; j<mu2e::CrvDigi::NSamples; j++) {
	      std::cout << cc.DBV_ADCs(pos,i)[j];
	      if(j<mu2e::CrvDigi::NSamples-1) {
		std::cout << " ";
	      }
	    }
	    std::cout << std::endl;
         

	    std::cout << "timestamp: " << timestamp << std::endl;
	    std::cout << "sysID: " << sysID << std::endl;
	    std::cout << "dtcID: " << dtcID << std::endl;
	    std::cout << "rocID: " << rocID << std::endl;
	    std::cout << "packetCount: " << packetCount << std::endl;
	    std::cout << "valid: " << valid << std::endl;
	    std::cout << "EVB mode: " << EVBMode << std::endl;
	    
//	    for(int i=7; i>=0; i--) {
//	      std::cout << (adc_t) *(pos+8+i);
//	      std::cout << " ";
//	    }
//	    std::cout << std::endl;
//	    
//	    for(int i=7; i>=0; i--) {
//	      std::cout << (adc_t) *(pos+8*2+i);
//	      std::cout << " ";
//	    }
//	    std::cout << std::endl;
	  
	    std::cout << "SiPMNumber: " << cc.DBV_sipmID(pos,i)%4 << std::endl;
	    std::cout << "scintillatorBarIndex: " << cc.DBV_sipmID(pos,i)/4 << std::endl;
	    std::cout << "TDC: " << cc.DBV_startTDC(pos,i) << std::endl;
	    std::cout << "Waveform: {";
	    for(size_t j=0; j<mu2e::CrvDigi::NSamples; j++) {
	      std::cout << cc.DBV_ADCs(pos,i)[j];
	      if(j<mu2e::CrvDigi::NSamples-1) {
		std::cout << ",";
	      }
	    }
	    std::cout << "}" << std::endl;
	  }

	  
	  std::cout << "LOOP: " << eventNumber << " " << curBlockIdx << " " << "(" << timestamp << ")" << std::endl;	    

	  for(size_t i=0; i<cc.DBVR_NumHits(pos); i++) {
	    // Text format: timestamp sipmID tdc nsamples sample_list
	    std::cout << "GREPMECRV: " << timestamp << " ";
	    std::cout << cc.DBV_sipmID(pos,i) << " ";
	    std::cout << cc.DBV_startTDC(pos,i) << " ";
	    for(size_t j=0; j<mu2e::CrvDigi::NSamples; j++) {
	      std::cout << cc.DBV_ADCs(pos,i)[j];
	      if(j<mu2e::CrvDigi::NSamples-1) {
		std::cout << " ";
	      }
	    }
	    std::cout << std::endl;
	  }

	} // End debug output
      } // End parsing CRV packets
    } // End loop over DataBlocks within fragment 
  } // Close loop over fragments

  if( diagLevel_ > 0 ) {
    std::cout << "mu2e::CrvDigisFromFragments::produce exiting eventNumber=" << (int)(event.event()) << " / timestamp=" << (int)eventNumber <<std::endl;
  }

  event.put(std::unique_ptr<EventNumber_t>(new EventNumber_t( eventNumber )));
  
  // Store the straw digis and calo digis in the event
  event.put(std::move(crv_digis));

}  // produce()

// ======================================================================

DEFINE_ART_MODULE(CrvDigisFromFragments)

// ======================================================================
