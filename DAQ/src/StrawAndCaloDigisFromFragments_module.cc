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
#include "mu2e-artdaq-core/Overlays/FragmentType.hh"
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
  struct Config 
  {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    fhicl::Atom<int>            diagLevel             { Name("diagLevel"),         Comment("diagnostic level")};
    fhicl::Atom<int>            parseCAL              { Name("parseCAL"),          Comment("parseCAL")};
    fhicl::Atom<int>            parseTRK              { Name("parseTRK"),          Comment("parseTRK")};
    fhicl::Atom<art::InputTag>  caloTag               { Name("caloTag"),           Comment("caloTag") };
    fhicl::Atom<art::InputTag>  trkTag                { Name("trkTag"),            Comment("trkTag") };
  };
  using EventNumber_t = art::EventNumber_t;
  using adc_t = mu2e::ArtFragmentReader::adc_t;
  
  // --- C'tor/d'tor:
  explicit  StrawAndCaloDigisFromFragments(const art::EDProducer::Table<Config>& config);
  virtual  ~StrawAndCaloDigisFromFragments()  { }

  // --- Production:
  virtual void produce( Event & );

private:
  int   diagLevel_;

  int   parseCAL_;
  int   parseTRK_;

  art::InputTag trkFragmentsTag_;
  art::InputTag caloFragmentsTag_;
  
  const int hexShiftPrint = 7;

};  // StrawAndCaloDigisFromFragments

// ======================================================================

StrawAndCaloDigisFromFragments::StrawAndCaloDigisFromFragments(const art::EDProducer::Table<Config>& config):
  art::EDProducer{ config },
  diagLevel_       (config().diagLevel()),
  parseCAL_        (config().parseCAL()),
  parseTRK_        (config().parseTRK()),
  trkFragmentsTag_ (config().trkTag()),
  caloFragmentsTag_(config().caloTag()){
    if (parseTRK_){
      produces<mu2e::StrawDigiCollection>();
    }
    if (parseCAL_){
      produces<mu2e::CaloDigiCollection>();
    }
  }

// ----------------------------------------------------------------------

void
StrawAndCaloDigisFromFragments::
produce( Event & event )
{
  art::EventNumber_t eventNumber = event.event();

  // Collection of StrawDigis for the event
  std::unique_ptr<mu2e::StrawDigiCollection> straw_digis(new mu2e::StrawDigiCollection);

  // Collection of CaloDigis for the event
  std::unique_ptr<mu2e::CaloDigiCollection> calo_digis(new mu2e::CaloDigiCollection);

  art::Handle<artdaq::Fragments> trkFragments, calFragments;
  size_t numTrkFrags(0), numCalFrags(0);
  if (parseTRK_){
    event.getByLabel(trkFragmentsTag_ , trkFragments);
    if (!trkFragments.isValid()){       
      std::cout << "[StrawAndCaloDigisFromFragments::produce] found no Tracker fragments!" << std::endl;
      event.put(std::move(straw_digis));
      return;
    }
    numTrkFrags = trkFragments->size();
  }
  if (parseCAL_){
    event.getByLabel(caloFragmentsTag_, calFragments);
    if (!calFragments.isValid()){
      std::cout << "[StrawAndCaloDigisFromFragments::produce] found no Calorimeter fragments!" << std::endl;
      event.put(std::move(calo_digis));
      return;
    }
    numCalFrags = calFragments->size();
  }
  // size_t numTrkFrags = trkFragments->size();
  // size_t numCalFrags = calFragments->size();

  if( diagLevel_ > 1 ) {
    std::cout << std::dec << "Producer: Run " << event.run() << ", subrun " << event.subRun()
	      << ", event " << eventNumber << " has " << std::endl;
    std::cout << numTrkFrags << " TRK fragments, and ";
    std::cout << numCalFrags << " CAL fragments." << std::endl;

    size_t totalSize = 0;
    for(size_t idx = 0; idx < numTrkFrags; ++idx) {
      auto size = ((*trkFragments)[idx]).size();// * sizeof(artdaq::RawDataType);
      totalSize += size;
      //      std::cout << "\tTRK Fragment " << idx << " has size " << size << std::endl;
    }
    for(size_t idx = 0; idx < numCalFrags; ++idx) {
      auto size = ((*calFragments)[idx]).size();// * sizeof(artdaq::RawDataType);
      totalSize += size;
      //      std::cout << "\tCAL Fragment " << idx << " has size " << size << std::endl;
    }
    totalSize *= sizeof(artdaq::RawDataType);

    std::cout << "\tTotal Size: " << (int)totalSize << " bytes." << std::endl;  
  }
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
    
    mu2e::FragmentType mode_;

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
	for(int i=hexShiftPrint; i>=0; i--) {
	  std::cout <<"0x" << std::hex << std::setw(4) << std::setfill('0')<< (adc_t) *(pos+i) << std::dec << std::setw(0);
	  std::cout << " ";
	}
	std::cout << std::endl;
      }	    

      auto hdr = cc.GetHeader(curBlockIdx);
      if(hdr == nullptr) {
	mf::LogError("StrawAndCaloDigisFromFragments") << "Unable to retrieve header from block " << curBlockIdx << "!" << std::endl;
	continue;
      }

      if(diagLevel_ > 1) {


	std::cout << "timestamp: " << static_cast<int>(hdr->GetTimestamp()) << std::endl;
	std::cout << "hdr->SubsystemID: " << static_cast<int>(hdr->SubsystemID) << std::endl;
	std::cout << "dtcID: " << static_cast<int>(hdr->DTCID) << std::endl;
	std::cout << "rocID: " << static_cast<int>(hdr->ROCID) << std::endl;
	std::cout << "packetCount: " <<static_cast<int>( hdr->PacketCount) << std::endl;
	std::cout << "valid: " << static_cast<int>(hdr->Valid) << std::endl;
	std::cout << "EVB mode: " << static_cast<int>(hdr->EVBMode) << std::endl;

	    for(int i=hexShiftPrint; i>=0; i--) {
	      std::cout << (adc_t) *(pos+8+i);
	      std::cout << " ";
	    }
	    std::cout << std::endl;
      }

      eventNumber = hdr->GetTimestamp();
      
      if(idx < numTrkFrags){
	mode_ = mu2e::FragmentType::TRK;//"TRK";
      }else {
	mode_ = mu2e::FragmentType::CAL;//"CAL";
      }

      // Parse phyiscs information from TRK packets
      if(mode_ == mu2e::FragmentType::TRK && hdr->PacketCount>0 && parseTRK_>0) {

	// Create the StrawDigi data products
	auto trkData = cc.GetTrackerData(curBlockIdx);
	if(trkData == nullptr) {
	  mf::LogError("StrawAndCaloDigisFromFragments") << "Error retrieving Tracker data from DataBlock " << curBlockIdx << "! Aborting processing of this block!";
	  continue;
	}

	mu2e::StrawId sid(trkData->StrawIndex);
	mu2e::TrkTypes::TDCValues tdc = {trkData->TDC0 , trkData->TDC1};
	mu2e::TrkTypes::TOTValues tot = {trkData->TOT0 , trkData->TOT1};


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


	mu2e::TrkTypes::ADCWaveform wf = trkData->Waveform();	

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
	  	  
	  for(int i=hexShiftPrint; i>=0; i--) {
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
	    if( ((0x0001<<(15-i)) & trkData->PreprocessingFlags) > 0) {
	      std::cout << "1";
	    } else {
	      std::cout << "0";
	    }
	  }
	  std::cout << std::endl;
	  	  
	  std::cout << "LOOP: " << eventNumber << " " << curBlockIdx << " " << "(" << hdr->GetTimestamp() << ")" << std::endl;	    
	  
	  // Text format: timestamp strawidx tdc0 tdc1 nsamples sample0-11
	  // Example: 1 1113 36978 36829 12 1423 1390 1411 1354 2373 2392 2342 2254 1909 1611 1525 1438
	  std::cout << "GREPMETRK: " << hdr->GetTimestamp() << " ";
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


      } else if(mode_ == mu2e::FragmentType::CAL && hdr->PacketCount>0 && parseCAL_>0) {	// Parse phyiscs information from CAL packets
	
	auto calData = cc.GetCalorimeterData(curBlockIdx);
	if(calData == nullptr) {
	  mf::LogError("StrawAndCaloDigisFromFragments") << "Error retrieving Calorimeter data from block " << curBlockIdx << "! Aborting processing of this block!";
	  continue;
	}

	if( diagLevel_ > 0 ) {
	  std::cout <<"[StrawAndCaloDigiFromFragments] NEW CALDATA: NumberOfHits "<< calData->NumberOfHits << std::endl;
	}

	bool err = false;
	for(size_t hitIdx = 0; hitIdx<calData->NumberOfHits; hitIdx++) {

	  // Fill the CaloDigiCollection
	  const mu2e::ArtFragmentReader::CalorimeterHitReadoutPacket* hitPkt(0);
	  hitPkt = cc.GetCalorimeterReadoutPacket(curBlockIdx, hitIdx);
	  if(hitPkt == nullptr) {
	    mf::LogError("StrawAndCaloDigisFromFragments") << "Error retrieving Calorimeter data from block " << curBlockIdx << " for hit " << hitIdx << "! Aborting processing of this block!";
	    err = true;
	    break;
	  }
	  
	  if( diagLevel_ > 0 ) {
	    std::cout <<"[StrawAndCaloDigiFromFragments] calo hit "<< hitIdx <<std::endl;
	    std::cout <<"[StrawAndCaloDigiFromFragments] \thitPkt " << hitPkt << std::endl;
	    std::cout <<"[StrawAndCaloDigiFromFragments] \tChNumber   " << (int)hitPkt->ChannelNumber  << std::endl;
	    std::cout <<"[StrawAndCaloDigiFromFragments] \tDIRACA     " << (int)hitPkt->DIRACA 	  << std::endl;
	    std::cout <<"[StrawAndCaloDigiFromFragments] \tDIRACB     " << (int)hitPkt->DIRACB	  << std::endl;
	    std::cout <<"[StrawAndCaloDigiFromFragments] \tErrorFlags " << (int)hitPkt->ErrorFlags	  << std::endl;
	    std::cout <<"[StrawAndCaloDigiFromFragments] \tTime	      " << (int)hitPkt->Time		  << std::endl;
	    std::cout <<"[StrawAndCaloDigiFromFragments] \tNSamples   " << (int)hitPkt->NumberOfSamples << std::endl;
	    std::cout <<"[StrawAndCaloDigiFromFragments] \tIndexMax   " << (int)hitPkt->IndexOfMaxDigitizerSample << std::endl;
	  }
	  
	  auto first = cc.GetCalorimeterReadoutSample(curBlockIdx,hitIdx,0);
	  auto last  = cc.GetCalorimeterReadoutSample(curBlockIdx, hitIdx, hitPkt->NumberOfSamples - 1);
	  if(first == nullptr || last == nullptr) {
	    mf::LogError("StrawAndCaloDigisFromFragments") << "Error retrieving Calorimeter samples from block " << curBlockIdx << " for hit " << hitIdx << "! Aborting processing of this block!";
	    err = true;
	    break;
	  }

	  //the second argument is not included in the vector, so we need to add "+1"
	  // because we want the "last" item included
	  std::vector<int> cwf(first,last+1);

	  // IMPORTANT NOTE: we don't have a final
	  // mapping yet so for the moment, the BoardID field (described in docdb 4914) is just a
	  // placeholder. Because we still need to know which crystal a hit belongs to, we are
	  // temporarily storing the 4-bit apdID and 12-bit crystalID in the Reserved DIRAC A slot.
	  // Also, note that until we have an actual map, channel index does not actually correspond
	  // to the physical readout channel on a ROC.
	  adc_t crystalID  = hitPkt->DIRACB & 0x0FFF;
	  adc_t apdID      = hitPkt->DIRACB >> 12;

	  calo_digis->emplace_back((crystalID*2 + apdID),
				   hitPkt->Time,
				   cwf,
				   hitPkt->IndexOfMaxDigitizerSample
				   );

	  if( diagLevel_ > 1 ) {
	    // Until we have the final mapping, the BoardID is just a placeholder
	    // adc_t BoardId    = cc.DBC_BoardID(pos,channelIdx);
        	      
	    std::cout << "Crystal ID: " << (int)crystalID << std::endl;		
	    std::cout << "APD ID: " << (int)apdID << std::endl;
	    std::cout << "Time: " << (int)hitPkt->Time << std::endl;
	    std::cout << "NumSamples: " << (int)hitPkt->NumberOfSamples << std::endl;
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
	    std::cout << "GREPMECAL: " << hdr->GetTimestamp() << " ";
	    std::cout << crystalID << " ";
	    std::cout << apdID << " ";
	    std::cout << hitPkt->Time << " ";
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
	if(err) continue;
	
      } // End Cal Mode
      
    } // End loop over DataBlocks within fragment 
      
  } // Close loop over fragments

  if( diagLevel_ > 0 ) {
    std::cout << "mu2e::StrawAndCaloDigisFromFragments::produce exiting eventNumber=" << (int)(event.event()) << " / timestamp=" << (int)eventNumber <<std::endl;

  }

  // Store the straw digis and calo digis in the event
  if (parseTRK_){
    event.put(std::move(straw_digis));
  }
  if (parseCAL_){
    event.put(std::move(calo_digis));
  }

}  // produce()

// ======================================================================

DEFINE_ART_MODULE(StrawAndCaloDigisFromFragments)

// ======================================================================
