//
// EDProducer module for converting calorimeter digis into DTC formatted
// packets (which are stored in DataBlock data products)
//
//
// Original author Tomonari Miyashita
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

// Mu2e includes.
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloDigiCollection.hh"
#include "DAQDataProducts/inc/DataBlockCollection.hh"

#include "SeedService/inc/SeedService.hh"

#include <fstream>
#include <stdexcept>

using namespace std;

namespace mu2e {

  //--------------------------------------------------------------------
  //
  //
  class CaloPacketProducer : public art::EDProducer {
  public:

    using adc_t = mu2e::DataBlock::adc_t;
    using dtc_id = mu2e::DataBlock::dtc_id;

    struct calhit {
      DataBlock::timestamp evt;
      size_t crystalId;
      int apdID;
      int recoDigiId;
      int recoDigiT0;
      int recoDigiSamples;
      std::vector<adc_t> waveform;

      dtc_id dtcID;
      int rocID;
    };

    explicit CaloPacketProducer(fhicl::ParameterSet const& pset);

    virtual void beginJob() override;

    virtual void endJob();

    void produce( art::Event & ) override;

  private:

    string                _outputFile;
    ofstream              outputStream;


    const size_t number_of_rocs = 172;
    const size_t number_of_crystals_per_roc = 8;
    const size_t number_of_rocs_per_dtc = 6;

    // 6 rocs per DTC => 27 DTCs
    // 172 rocs * 8 crystals per roc => 1376
    // Note: the highest crystal ID in the old simulation was 1355

    int _generateTextFile;

    // Diagnostics level.
    int _diagLevel;

    // Enable DIRAC emulation
    int _enableDIRACEmulation;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    // Label of the module that made the digis.
    std::string _makerModuleLabel;

  };

  CaloPacketProducer::CaloPacketProducer(fhicl::ParameterSet const& pset):
    art::EDProducer{pset},
    _outputFile                    (pset.get<string>("outputFile","artdaq_cal.txt")),
    _generateTextFile(pset.get<int>("generateTextFile",0)),
    _diagLevel(pset.get<int>("diagLevel",0)),
    _enableDIRACEmulation(pset.get<int>("enableDIRACEmulation",0)),
    _maxFullPrint(pset.get<int>("maxFullPrint",5)),
    _makerModuleLabel(pset.get<std::string>("makerModuleLabel","CaloDigiFromShower")) {

    produces<DataBlockCollection>();

    if(_generateTextFile>0) {
      outputStream.open(_outputFile, std::ofstream::out | std::ofstream::trunc);
    }
  }

  void CaloPacketProducer::beginJob(){

    if ( _diagLevel > 0 ) {
      cout << "CaloPacketProducer Diaglevel: "
           << _diagLevel << " "
           << _maxFullPrint
           << endl;
    }

  }

  void CaloPacketProducer::endJob(){
    if(_generateTextFile>0) {
      outputStream << flush;
      outputStream.close();
    }
  }


  void CaloPacketProducer::produce(art::Event & evt) {

    unique_ptr<DataBlockCollection> dtcPackets(new DataBlockCollection);

    bool gblresult;
    art::Handle<CaloDigiCollection> caloDigisHandle;
    gblresult = evt.getByLabel(_makerModuleLabel, caloDigisHandle);
    ( _diagLevel > 0 ) &&
      std::cout
      << __func__ << " getting data by getByLabel: label, instance, result " << std::endl
      << " CaloDigiCollection             _makerModuleLabel: " << _makerModuleLabel << ", "
      << ", "<< gblresult << std::endl;

    if (!gblresult) throw cet::exception("DATA") << " Missing digi data";

    CaloDigiCollection const& hits_CD = gblresult ? *caloDigisHandle : CaloDigiCollection();

    auto eventNum = evt.id().event();
    if ( _diagLevel > 2 ) {
      cout << "CaloPacketProducer: eventNum: " << eventNum
	   << " Total number of calo hit digis = " << hits_CD.size() << endl;
    }

    std::vector<calhit> caloHitVector; // Vector of calo hit digi data
    for ( size_t i=0; i<hits_CD.size(); ++i ) {

      CaloDigi const& CD = gblresult ? hits_CD.at(i) : CaloDigi();
      const std::vector<int> theWaveform = CD.waveform();

      // Fill struct with info for current hit
      calhit curHit;
      curHit.evt = eventNum;
      curHit.crystalId = CD.roId()/2;
      curHit.recoDigiId = CD.roId();
      curHit.recoDigiT0 = CD.t0();
      curHit.recoDigiSamples = theWaveform.size();
      for(size_t j = 0; j<theWaveform.size(); j++) {
	curHit.waveform.push_back(theWaveform[j]);
      }

      // 192 ROCs total
      if(curHit.crystalId >= number_of_rocs * number_of_crystals_per_roc) {
	throw cet::exception("DATA") << " Crystal index " << curHit.crystalId
				     << " exceeds limit of " << number_of_rocs
				     << " * " << number_of_crystals_per_roc
				     << " =" << number_of_rocs*number_of_crystals_per_roc;
      }

      // ROC ID, counting from 0, across all (for the calorimeter)
      size_t globalROCID = curHit.crystalId / number_of_crystals_per_roc;

      curHit.rocID = globalROCID % number_of_rocs_per_dtc;
      curHit.apdID = curHit.recoDigiId % 2; // Even is APD 0, Odd is APD 1
      curHit.dtcID = dtc_id(globalROCID / number_of_rocs_per_dtc);

      caloHitVector.push_back(curHit);

      if(_generateTextFile>0) {
	outputStream << curHit.evt << "\t";
	outputStream << curHit.crystalId << "\t";
	//	outputStream << curHit.recoDigiId << "\t"; // Readout ID
	outputStream << curHit.apdID << "\t"; // Readout ID
	outputStream << curHit.recoDigiT0 << "\t"; // Readout time
	outputStream << curHit.recoDigiSamples << "\t";
	for(size_t j=0; j<curHit.waveform.size(); j++) {
	  outputStream << curHit.waveform[j];
	  if(j<curHit.waveform.size()-1) {
	    outputStream << "\t";
	  }
	}
	outputStream << endl;
      }
      
    }


    dtc_id max_dtc_id = number_of_rocs/number_of_rocs_per_dtc-1;
    if(number_of_rocs % number_of_rocs_per_dtc > 0) {
      max_dtc_id += 1;
    }

    // Create a vector to hold all DTC ID / ROC ID combinations
    std::vector< std::pair<dtc_id, size_t> > dtcRocVector;
    for(dtc_id curDTCID = 0; curDTCID <= max_dtc_id; curDTCID++) {
      for(size_t curROCID = 0; curROCID < number_of_rocs_per_dtc; curROCID++) {
	std::pair<dtc_id, size_t> curPair(curDTCID, curROCID);
	dtcRocVector.push_back(curPair);	  
      }
    }

    // // Randomize the order in which the DTCs are added to the event
    // std::shuffle(dtcRocVector.begin(),dtcRocVector.end(),generator);

    // Loop over the DTC/ROC pairs and generate datablocks for each ROC
    auto targetNumROCs = dtcRocVector.size();
    for(size_t curPairNum = 0; curPairNum < targetNumROCs; curPairNum++) {

      size_t dtcID = dtcRocVector[curPairNum].first;
      size_t rocID = dtcRocVector[curPairNum].second;

      std::vector<calhit> curHitVector;
      // Find all hits for this event coming from the specified DTC/ROC
      for (size_t curHitIdx = 0; curHitIdx < caloHitVector.size(); curHitIdx++) {
	if (caloHitVector[curHitIdx].dtcID == (int)dtcID &&
	    caloHitVector[curHitIdx].rocID == (int)rocID) {
	  curHitVector.push_back(caloHitVector[curHitIdx]);
	}
      }

      if (curHitVector.size() == 0) {

	// No hits
	std::vector<adc_t> curDataBlock;
	adc_t null_adc = 0;
	// First 16 bits of header (reserved values)
	curDataBlock.push_back(null_adc);
	// Second 16 bits of header (ROC ID, packet type):
	adc_t curROCID = rocID; // 4 bit ROC ID
	adc_t headerPacketType = 5; // 4 bit Data packet header type is 5
	headerPacketType <<= 4; // Shift left by 4
	adc_t secondEntry = (curROCID | headerPacketType);
	secondEntry = (secondEntry | (1 << 15)); // valid bit
	curDataBlock.push_back(secondEntry);
	// Third 16 bits of header (the number of data packets following the header
	// is 1 when there are no hits)
	curDataBlock.push_back(0x0001);
	// Fourth through sixth 16 bits of header (timestamp)
	uint64_t timestamp = eventNum;
	curDataBlock.push_back(static_cast<adc_t>(timestamp & 0xFFFF));
	curDataBlock.push_back(static_cast<adc_t>((timestamp >> 16) & 0xFFFF));
	curDataBlock.push_back(static_cast<adc_t>((timestamp >> 32) & 0xFFFF));
	
	// Seventh 16 bits of header (data packet format version and status)
	adc_t status = 0; // 0 Corresponds to "Timestamp has valid data"
	adc_t formatVersion = (5 << 8); // Using 5 for now
	curDataBlock.push_back(formatVersion + status);

	// Eighth 16 bits of header (EVB Mode | SYSID | DTCID)
	adc_t evbMode = (0 << 8);
	adc_t sysID = (1 << 6) & 0x00C0;
	adc_t curDTCID = dtcID & 0x003F;
	curDataBlock.push_back(evbMode + sysID + curDTCID);
	
	// Fill in the byte count field of the header packet
	adc_t numBytes = 16*2;
	// A ROC with no hits requires 2 packets (including the header packet)
	curDataBlock[0] = numBytes;

	///////////////////////////////////////////
	// Create the CAL ROC Hit Index Packet
	///////////////////////////////////////////
	adc_t numHits = 0x0000;
	curDataBlock.push_back(numHits);

	// Add the Board ID and Channel Status fields. For now, we'll just fill this with 0x0000.
	// Once accessor methods and a final conversion table are available, this should be
	// replaced with an actual board ID as described in docdb 4914. For now, we'll also set
	// all the Channel Status Flags to 1. Once the full wiring map is available, this can be
	// updated to take into account cases where fewer than 20 channels are connected to a ROC.
	curDataBlock.push_back(0xFC00);
	curDataBlock.push_back(0x3FFF);

	// Pad the empty space in the last packet. 3 fields were used for the number of hits,
	// the Board ID, and the Channel Status Flags so there are 5 fields left in the current
	// packet.
	for(size_t i=0; i<5; i++) {
	  curDataBlock.push_back(0x0000);
	}
	
	// Make sure we end on a packet boundary
	assert(curDataBlock.size() % 8 == 0);

	/////////////////////////////////////////////////////
	// Create mu2e::DataBlock and add to the collection
	/////////////////////////////////////////////////////	
	DataBlock theBlock(DataBlock::CAL, evt.id(), dtcID, curDataBlock);
	dtcPackets->push_back(theBlock);
	
      } else {

	///////////////////////////////////////////
	// Create the header packet
	///////////////////////////////////////////

	std::vector<adc_t> curDataBlock;
	// Add the header packet to the DataBlock (leaving including a placeholder for
	// the number of packets in the DataBlock);
	adc_t null_adc = 0;
	// First 16 bits of header (reserved values)
	curDataBlock.push_back(null_adc);
	// Second 16 bits of header (ROC ID, packet type):
	adc_t curROCID = rocID; // 4 bit ROC ID
	adc_t headerPacketType = 5; // 4 bit Data packet header type is 5
	headerPacketType <<= 4; // Shift left by 4
	adc_t secondEntry = (curROCID | headerPacketType);
	secondEntry = (secondEntry | (1 << 15)); // valid bit
	curDataBlock.push_back(secondEntry);
	// Third 16 bits of header (number of data packets is 0)
	curDataBlock.push_back(null_adc);
	// Fourth through sixth 16 bits of header (timestamp)
	uint64_t timestamp = eventNum;
	curDataBlock.push_back(static_cast<adc_t>(timestamp & 0xFFFF));
	curDataBlock.push_back(static_cast<adc_t>((timestamp >> 16) & 0xFFFF));
	curDataBlock.push_back(static_cast<adc_t>((timestamp >> 32) & 0xFFFF));
	
	// Seventh 16 bits of header (data packet format version and status)
	adc_t status = 0; // 0 Corresponds to "Timestamp has valid data"
	adc_t formatVersion = (5 << 8); // Using 5 for now
	curDataBlock.push_back(formatVersion + status);

	// Eighth 16 bits of header (EVB Mode | SYSID | DTCID)
	adc_t evbMode = (0 << 8);
	adc_t sysID = (1 << 6) & 0x00C0;
	adc_t curDTCID = dtcID & 0x003F;
	curDataBlock.push_back(evbMode + sysID + curDTCID);
	
	///////////////////////////////////////////
	// Create the CAL ROC Hit Index Packets
	///////////////////////////////////////////

	size_t numHits = curHitVector.size();
	
	// The first field in the CAL ROC Hit Index Packet is the number of hits
	curDataBlock.push_back((adc_t)numHits);

	// The first offset is always numHits+3 because there is 1 field for the number
	// of hits, numHits fields for the hit indices, and 2 fields for the board ID
	// and Channel Status.
	size_t curOffset = numHits+3;
	curDataBlock.push_back(curOffset);
	for(size_t curHitIdx = 1; curHitIdx < numHits; curHitIdx++) {
	  calhit previousHit = curHitVector[curHitIdx-1];
	  size_t previousHitSize = previousHit.waveform.size() + 5;
	  // In addition to the waveform samples, there are 5 other fields in a channel readout:
	  // DIRAC A, DIRAC B, Error Flags, Time, and the Num Samples/Max Sample Idx field
	  curOffset += previousHitSize;
	  curDataBlock.push_back(curOffset);
	}
	
	// Find the offset of the end of the last hit
	curOffset += curHitVector[numHits-1].waveform.size() + 5;

	// Add the Board ID and Channel Status fields. For now, we'll just fill this with 0x0000.
	// Once accessor methods and a final conversion table are available, this should be
	// replaced with an actual board ID as described in docdb 4914. For now, we'll also set
	// all the Channel Status Flags to 1. Once the full wiring map is available, this can be
	// updated to take into account cases where fewer than 20 channels are connected to a ROC.
	curDataBlock.push_back(0xFC00);
	curDataBlock.push_back(0x3FFF);

	///////////////////////////////////////////
	// Create the channel hit list
	///////////////////////////////////////////

	for(size_t curHitIdx = 0; curHitIdx<numHits; curHitIdx++) {
	  calhit curHit = curHitVector[curHitIdx];

	  // Assume the 0th apd is always read out before the second
	  adc_t crystalID = curHit.crystalId;
	  adc_t apdID = curHit.apdID;
	  adc_t IDNum = ((apdID << 12) | crystalID);

	  // Fill the channel number and Reserved For DIRAC A fields
	  // IMPORTANT NOTE: Since we don't have the wiring map, we'll just fill this with 0s for now
	  // Also, the Reserved For DIRAC A field is not yet defined, so we'll fill it with 0s as well
	  adc_t reservedForDIRACA = 0x0000;
	  adc_t channelNum = 0x0000;
	  adc_t firstField = (reservedForDIRACA & 0xFFC0) | (channelNum & 0x003F);
	  curDataBlock.push_back(firstField);

	  // Fill field second field reserved for DIRAC output
	  // IMPORTANT NOTE: Since we are not using actual Board IDs (see comment above) we still need
	  // a way to store and recover the crystal ID until we have an actual wiring map. For now,
	  // we will simply store the apdID and crystalID in the first per-hit field reserved for DIRAC
	  // output
	  curDataBlock.push_back(IDNum);

	  // Fill Error Flag field (format TBD):
	  curDataBlock.push_back(0x0000);

	  // Fill Time field:
	  curDataBlock.push_back((adc_t)(curHit.recoDigiT0));
	  
	  // Fill Num Samples / DIRAC Max Peak Idx Field:
	  if(_enableDIRACEmulation==0) {
	    curDataBlock.push_back((adc_t)(curHit.recoDigiSamples));
	  } else {
	    adc_t maxVal = 0;
	    size_t maxIdx = 0;
	    for (auto sampleIdx = 0; sampleIdx < curHit.recoDigiSamples; sampleIdx++) {
	      adc_t scaledVal = static_cast<adc_t>(curHit.waveform[sampleIdx]);
	      if(scaledVal>=maxVal) {
		maxVal = scaledVal;
		maxIdx = sampleIdx;
	      }
	    }
	    curDataBlock.push_back( (0xFF00 & (adc_t)(maxIdx<<8)) | (0x00FF & (adc_t)(curHit.recoDigiSamples)) );
	  }
	  
	  for (auto sampleIdx = 0; sampleIdx < curHit.recoDigiSamples; sampleIdx++) {
	      adc_t scaledVal = static_cast<adc_t>(curHit.waveform[sampleIdx]);
	      curDataBlock.push_back(scaledVal);
	  }

	} // Done looping over hits for this DTC/ROC pair


//	for(size_t curHitIdx=0; curHitIdx<numHits; curHitIdx++) {
//	  calhit curHit = curHitVector[curHitIdx];
//	  std::cout << "\t" << dtcID << "/" << rocID << "\t" << curHit.crystalId << "\t" << std::dec << curHit.dtcID << "/" << curHit.rocID << std::endl;
//	  std::cout << std::endl;
//	}
	
	// Pad any empty space in the last packet with 0s
	size_t padding_slots = 8 - (curDataBlock.size() % 8);
	if(padding_slots<8) {
	    for (size_t i = 0; i < padding_slots; i++) {
	      curDataBlock.push_back(0x0000);
	    }
	}

	// Fill in the number of data packets entry in the header packet
	adc_t numDataPackets = static_cast<adc_t>(curDataBlock.size() / 8 - 1);
	curDataBlock[2] = numDataPackets;

	// Fill in the byte count field of the header packet
	adc_t numBytes = curDataBlock.size() * 2;
	curDataBlock[0] = numBytes;

	// Make sure we end on a packet boundary
	assert(curDataBlock.size() % 8 == 0);
	

	/////////////////////////////////////////////////////
	// Create mu2e::DataBlock and add to the collection
	/////////////////////////////////////////////////////	
	DataBlock theBlock(DataBlock::CAL, evt.id(), dtcID, curDataBlock);
	dtcPackets->push_back(theBlock);
      }

    } // Done looping over  DTC/ROC pairs
 
    evt.put(move(dtcPackets));

  } // end of ::produce

}


using mu2e::CaloPacketProducer;
DEFINE_ART_MODULE(CaloPacketProducer);
