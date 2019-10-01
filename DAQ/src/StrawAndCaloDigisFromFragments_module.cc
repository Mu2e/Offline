// ======================================================================
//
// StrawAndCaloDigisFromFragments_plugin:  Add tracker/cal data products to the event
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
#include "DataProducts/inc/TrkTypes.hh"
#include "RecoDataProducts/inc/StrawDigiCollection.hh"
#include "RecoDataProducts/inc/CaloDigiCollection.hh"

#include <iostream>

#include <string>

#include <memory>

namespace art {
  class StrawAndCaloDigisFromFragments;
}

using art::StrawAndCaloDigisFromFragments;

// ======================================================================

class art::StrawAndCaloDigisFromFragments
  : public EDProducer
{

public:

  using EventNumber_t = art::EventNumber_t;
  using adc_t = mu2e::ArtFragmentReader::adc_t;
  
  // --- C'tor/d'tor:
  explicit  StrawAndCaloDigisFromFragments(fhicl::ParameterSet const& pset);
  virtual  ~StrawAndCaloDigisFromFragments()  { }

  // --- Production:
  virtual void produce( Event & );

private:
  int   diagLevel_;

  int   parseCAL_;
  int   parseTRK_;

  art::InputTag trkFragmentsTag_;
  art::InputTag caloFragmentsTag_;

};  // StrawAndCaloDigisFromFragments

// ======================================================================

StrawAndCaloDigisFromFragments::StrawAndCaloDigisFromFragments(fhicl::ParameterSet const& pset)
  : EDProducer(pset)
  , diagLevel_(pset.get<int>("diagLevel",0))
  , parseCAL_(pset.get<int>("parseCAL",1))
  , parseTRK_(pset.get<int>("parseTRK",1))
  , trkFragmentsTag_(pset.get<art::InputTag>("trkTag","daq:trk"))
  , caloFragmentsTag_(pset.get<art::InputTag>("caloTag","daq:calo"))
{
  produces<EventNumber_t>(); 
  produces<mu2e::StrawDigiCollection>();
  produces<mu2e::CaloDigiCollection>();
}

// ----------------------------------------------------------------------

void
  StrawAndCaloDigisFromFragments::
  produce( Event & event )
{

  art::EventNumber_t eventNumber = event.event();

  auto trkFragments = event.getValidHandle<artdaq::Fragments>(trkFragmentsTag_);
  auto calFragments = event.getValidHandle<artdaq::Fragments>(caloFragmentsTag_);
  size_t numTrkFrags = trkFragments->size();
  size_t numCalFrags = calFragments->size();

  if( diagLevel_ > 1 ) {
    std::cout << std::dec << "Producer: Run " << event.run() << ", subrun " << event.subRun()
	      << ", event " << eventNumber << " has " << std::endl;
    std::cout << trkFragments->size() << " TRK fragments, and ";
    std::cout << calFragments->size() << " CAL fragments." << std::endl;

    size_t totalSize = 0;
    for(size_t idx = 0; idx < trkFragments->size(); ++idx) {
      auto size = ((*trkFragments)[idx]).size() * sizeof(artdaq::RawDataType);
      totalSize += size;
      //      std::cout << "\tTRK Fragment " << idx << " has size " << size << std::endl;
    }
    for(size_t idx = 0; idx < calFragments->size(); ++idx) {
      auto size = ((*calFragments)[idx]).size() * sizeof(artdaq::RawDataType);
      totalSize += size;
      //      std::cout << "\tCAL Fragment " << idx << " has size " << size << std::endl;
    }
    
    std::cout << "\tTotal Size: " << (int)totalSize << " bytes." << std::endl;  
  }

  // Collection of StrawDigis for the event
  std::unique_ptr<mu2e::StrawDigiCollection> straw_digis(new mu2e::StrawDigiCollection);

  // Collection of CaloDigis for the event
  std::unique_ptr<mu2e::CaloDigiCollection> calo_digis(new mu2e::CaloDigiCollection);

  std::string curMode = "TRK";

  // Loop over the TRK and CAL fragments
  for (size_t idx = 0; idx < numTrkFrags+numCalFrags; ++idx) {

    auto curHandle = trkFragments;
    size_t curIdx = idx;
    if(idx>=numTrkFrags) {
      curIdx = idx-numTrkFrags;
      curHandle = calFragments;
    }
    const auto& fragment((*curHandle)[curIdx]);
    
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
    
    std::string mode_;

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
      
      if(sysID==0) {
	mode_ = "TRK";
      } else if(sysID==1) {
	mode_ = "CAL";
      }

      // Parse phyiscs information from TRK packets
      if(mode_ == "TRK" && packetCount>0 && parseTRK_>0) {

	// Create the StrawDigi data products
	mu2e::StrawId sid(cc.DBT_StrawIndex(pos));
	mu2e::TrkTypes::TDCValues tdc = {cc.DBT_TDC0(pos) , cc.DBT_TDC1(pos)};
	mu2e::TrkTypes::TOTValues tot = {cc.DBT_TOT0(pos) , cc.DBT_TOT1(pos)};


//	///////////////////////////////////////////////////////////////////////////
//	// NOTE: Because the tracker code in offline has not been updated to
//	// use 15 samples, it is necessary to add an extra sample in order to
//	// initialize an ADCWaveform that can be passed to the StrawDigi
//	// constructor. This means that the digis produced by StrawAndCaloDigisFromFragments
//	// will differ from those processed in offline so the filter performance
//	// will be different. This is only temporary.
//	std::array<adc_t,15> const & shortWaveform = cc.DBT_Waveform(pos);
//	mu2e::TrkTypes::ADCWaveform wf;
//	for(size_t i=0; i<15; i++) {
//	  wf[i] = shortWaveform[i];
//	}
//	wf[15] = 0;
//	///////////////////////////////////////////////////////////////////////////


	mu2e::TrkTypes::ADCWaveform wf = cc.DBT_Waveform(pos);	

	adc_t flags = cc.DBT_Flags(pos);

	// Fill the StrawDigiCollection
	straw_digis->emplace_back(sid, tdc, tot, wf);

	if( diagLevel_ > 1 ) {
  	  std::cout << "MAKEDIGI: " << sid.asUint16() << " " << tdc[0] << " " << tdc[1] << " "
	    << tot[0] << " " << tot[1] << " ";
	  for(size_t i=0; i<mu2e::TrkTypes::NADC; i++) {
	    std::cout << wf[i];
	    if(i<mu2e::TrkTypes::NADC-1) {
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
	  
	  for(int i=7; i>=0; i--) {
	    std::cout << (adc_t) *(pos+8+i);
	    std::cout << " ";
	  }
	  std::cout << std::endl;
	  
	  for(int i=7; i>=0; i--) {
	    std::cout << (adc_t) *(pos+8*2+i);
	    std::cout << " ";
	  }
	  std::cout << std::endl;
	  
	  std::cout << "strawIdx: " << sid.asUint16() << std::endl;
	  std::cout << "TDC0: " << tdc[0] << std::endl;
	  std::cout << "TDC1: " << tdc[1] << std::endl;
	  std::cout << "TOT0: " << tot[0] << std::endl;
	  std::cout << "TOT1: " << tot[1] << std::endl;
	  std::cout << "Waveform: {";
	  for(size_t i=0; i<mu2e::TrkTypes::NADC; i++) {
	    std::cout << wf[i];
	    if(i<mu2e::TrkTypes::NADC-1) {
	      std::cout << ",";
	    }
	  }
	  std::cout << "}" << std::endl;

	  std::cout << "FPGA Flags: ";
	  for(size_t i=8; i<16; i++) {
	    if( ((0x0001<<(15-i)) & flags) > 0) {
	      std::cout << "1";
	    } else {
	      std::cout << "0";
	    }
	  }
	  std::cout << std::endl;
	  	  
	  std::cout << "LOOP: " << eventNumber << " " << curBlockIdx << " " << "(" << timestamp << ")" << std::endl;	    
	  
	  // Text format: timestamp strawidx tdc0 tdc1 nsamples sample0-11
	  // Example: 1 1113 36978 36829 12 1423 1390 1411 1354 2373 2392 2342 2254 1909 1611 1525 1438
	  std::cout << "GREPMETRK: " << timestamp << " ";
	  std::cout << sid.asUint16() << " ";
	  std::cout << tdc[0] << " ";
	  std::cout << tdc[1] << " ";
	  std::cout << tot[0] << " ";
	  std::cout << tot[1] << " ";
	  std::cout << wf.size() << " ";
	  for(size_t i=0; i<mu2e::TrkTypes::NADC; i++) {
	    std::cout << wf[i];
	    if(i<mu2e::TrkTypes::NADC-1) {
	      std::cout << " ";
	    }
	  }
	  std::cout << std::endl;
	} // End debug output


      } else if(mode_ == "CAL" && packetCount>0 && parseCAL_>0) {	// Parse phyiscs information from CAL packets
	
	for(size_t hitIdx = 0; hitIdx<cc.DBC_NumHits(pos); hitIdx++) {

	  // Fill the CaloDigiCollection
	  std::vector<int> cwf = cc.DBC_Waveform(pos,hitIdx);

	  // IMPORTANT NOTE: As described in CaloPacketProducer_module.cc, we don't have a final
	  // mapping yet so for the moment, the BoardID field (described in docdb 4914) is just a
	  // placeholder. Because we still need to know which crystal a hit belongs to, we are
	  // temporarily storing the 4-bit apdID and 12-bit crystalID in the Reserved DIRAC A slot.
	  // Also, note that until we have an actual map, channel index does not actually correspond
	  // to the physical readout channel on a ROC.
	  calo_digis->emplace_back(
				   ((cc.DBC_DIRACOutputB(pos,hitIdx) & 0x0FFF)*2 +
				    (cc.DBC_DIRACOutputB(pos,hitIdx)>>12)),
				     cc.DBC_Time(pos,hitIdx),
				     cwf
				   );

	    if( diagLevel_ > 1 ) {
	      // Until we have the final mapping, the BoardID is just a placeholder
	      // adc_t BoardId    = cc.DBC_BoardID(pos,channelIdx);
	      
	      // As noted above, until we have the final mapping, the crystal ID and apdID
	      // will be temporarily stored in the Reserved DIRAC A slot.
	      adc_t crystalID  = cc.DBC_DIRACOutputB(pos,hitIdx) & 0x0FFF;
	      adc_t apdID      = cc.DBC_DIRACOutputB(pos,hitIdx) >> 12;
	      adc_t time       = cc.DBC_Time(pos,hitIdx);
	      adc_t numSamples = cc.DBC_NumSamples(pos,hitIdx);
	      //adc_t peakIdx    = cc.DBC_PeakSampleIdx(pos,hitIdx);
         
	      std::cout << "timestamp: " << timestamp << std::endl;
	      std::cout << "sysID: " << sysID << std::endl;
	      std::cout << "dtcID: " << dtcID << std::endl;
	      std::cout << "rocID: " << rocID << std::endl;
	      std::cout << "packetCount: " << packetCount << std::endl;
	      std::cout << "valid: " << valid << std::endl;
	      std::cout << "EVB mode: " << EVBMode << std::endl;		
         	  
	      for(int i=7; i>=0; i--) {
		std::cout << (adc_t) *(pos+8+i);
		std::cout << " ";
	      }
	      std::cout << std::endl;
         	  
	      for(int i=7; i>=0; i--) {
		std::cout << (adc_t) *(pos+8*2+i);
		std::cout << " ";
	      }
	      std::cout << std::endl;
	      
	      std::cout << "Crystal ID: " << crystalID << std::endl;		
	      std::cout << "APD ID: " << apdID << std::endl;
	      std::cout << "Time: " << time << std::endl;
	      std::cout << "NumSamples: " << numSamples << std::endl;
	      std::cout << "Waveform: {";
	      for(size_t i=0; i<cwf.size(); i++) {
		std::cout << cwf[i];
		if(i<cwf.size()-1) {
		  std::cout << ",";
		}
	      }
	      std::cout << "}" << std::endl;
         	  
	      // Text format: timestamp crystalID roID time nsamples samples...
	      // Example: 1 201 402 660 18 0 0 0 0 1 17 51 81 91 83 68 60 58 52 42 33 23 16
	      std::cout << "GREPMECAL: " << timestamp << " ";
	      std::cout << crystalID << " ";
	      std::cout << apdID << " ";
	      std::cout << time << " ";
	      std::cout << cwf.size() << " ";
	      for(size_t i=0; i<cwf.size(); i++) {
		std::cout << cwf[i];
		if(i<cwf.size()-1) {
		  std::cout << " ";
		}
	      }
	      std::cout << std::endl;
	    } // End debug output
	    
	} // End loop over readout channels in DataBlock
	
      } // End Cal Mode
      
    } // End loop over DataBlocks within fragment 
      
  } // Close loop over fragments

  //  }  // Close loop over the TRK and CAL collections

  if( diagLevel_ > 0 ) {
    std::cout << "mu2e::StrawAndCaloDigisFromFragments::produce exiting eventNumber=" << (int)(event.event()) << " / timestamp=" << (int)eventNumber <<std::endl;

  }

  event.put(std::unique_ptr<EventNumber_t>(new EventNumber_t( eventNumber )));
  
  // Store the straw digis and calo digis in the event
  event.put(std::move(straw_digis));
  event.put(std::move(calo_digis));

}  // produce()

// ======================================================================

DEFINE_ART_MODULE(StrawAndCaloDigisFromFragments)

// ======================================================================
