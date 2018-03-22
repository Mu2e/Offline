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
#include <math.h>

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
      int crystalId;
      int recoDigiId;
      int recoDigiT0;
      int recoDigiSamples;
      std::vector<adc_t> waveform;
      
      int rocID;
      int ringID;
      int apdID;
      
      dtc_id dtcID;
    };

    explicit CaloPacketProducer(fhicl::ParameterSet const& pset);

    virtual void beginJob() override;

    virtual void endJob();

    void produce( art::Event & ) override;

  private:

    string                _outputFile;
    ofstream              outputStream;

    // For now, crystal IDs start at 0 and end at number_of_rings*rocs_per_ring*number_of_crystals_per_roc-1
    const size_t number_of_crystals_per_roc = 8;

    // 6 optical transceivers per DTC
    // 2 optical links per ring (readout both ends)
    // => 3 Rings per DTC/server
    //
    // 16 ROCs/server 
    // 3 Rings per server
    // => 6 ROCs per ring (every third ring will actually only have 4 ROCS)
    //
    // 12 Servers
    // 3 rings per server
    // => 36 rings
    //
    // 12 servers
    // 16 ROCs/server
    // 8 crystals per ROC
    // => MAX 1536 Crystals
    //
    // However, for now if we use 6 rocs on every ring, we have 18 ROCs per
    // server giving us 216 ROCs for a maximum of 1728 crystals
    //
    const size_t number_of_rings = 36;
    const size_t rings_per_dtc = 3;
    const size_t rocs_per_ring = 6;
    // Actually, every third ring will only have 4 ROCs, but for now we'll
    // assign 6 to all rings

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
      curHit.crystalId = (int)(CD.roId()/2);
      curHit.recoDigiId = CD.roId();
      curHit.recoDigiT0 = CD.t0();
      curHit.recoDigiSamples = theWaveform.size();
      for(size_t j = 0; j<theWaveform.size(); j++) {
	curHit.waveform.push_back(theWaveform[j]);
      }

      // 192 ROCs total
      if(curHit.crystalId >= 192 * abs(number_of_crystals_per_roc)) {
	throw cet::exception("DATA") << " Crystal index " << curHit.crystalId
				     << " exceeds limit of 192*" << number_of_crystals_per_roc
				     << "=" << 216*number_of_crystals_per_roc;
      }

      // Ring ID, counting from 0, across all (for the calorimeter)
      size_t globalRingID = int(curHit.crystalId / (rocs_per_ring * number_of_crystals_per_roc));

      curHit.ringID = globalRingID % rings_per_dtc;
      curHit.rocID = (curHit.crystalId - (rocs_per_ring * number_of_crystals_per_roc) * globalRingID) / number_of_crystals_per_roc;
      curHit.apdID = curHit.recoDigiId % 2; // Even is APD 0, Odd is APD 1
      curHit.dtcID = dtc_id(globalRingID/rings_per_dtc);

      caloHitVector.push_back(curHit);

      if(_generateTextFile>0) {
	outputStream << curHit.evt << "\t";
	outputStream << curHit.crystalId << "\t";
	outputStream << curHit.recoDigiId << "\t"; // Readout ID
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
    
    dtc_id max_dtc_id = number_of_rings / rings_per_dtc - 1;
    if(number_of_rings % rings_per_dtc > 0) {
      max_dtc_id += 1;
    }
    // Create a vector to hold all Ring ID / ROC ID combinations for all DTCs to simulate
    std::vector< std::pair<dtc_id, std::pair<size_t,size_t> > > rocRingVector;
    for(dtc_id curDTCID = 0; curDTCID <= max_dtc_id; curDTCID++) {
      for(size_t curRingID = 0; curRingID < rings_per_dtc; curRingID++) {
	for(size_t curROCID = 0; curROCID < rocs_per_ring; curROCID++) {
	  std::pair<int, int> curRocRingPair(curRingID, curROCID);	  
	  std::pair<dtc_id, std::pair<size_t,size_t> > curPair(curDTCID, curRocRingPair);
	  rocRingVector.push_back(curPair);
	  
	}
      }
    }
    // // Randomize the order in which the Rings and ROCs are received
    // std::shuffle(rocRingVector.begin(),rocRingVector.end(),generator);


    // Loop over the ROC/ring pairs and generate datablocks for each ROC on
    // all the rings
    auto targetNumROCs = rocRingVector.size();
    for(size_t curPairNum = 0; curPairNum < targetNumROCs; curPairNum++) {

      size_t dtcID = rocRingVector[curPairNum].first;
      size_t ringID = rocRingVector[curPairNum].second.first;
      size_t rocID = rocRingVector[curPairNum].second.second;

      std::vector<calhit> curHitVector;
      // Find all hits for this event coming from the specified Ring/ROC
      for (size_t curHitIdx = 0; curHitIdx < caloHitVector.size(); curHitIdx++) {
	if (caloHitVector[curHitIdx].dtcID == (int)dtcID &&
	    caloHitVector[curHitIdx].rocID == (int)rocID &&
	    caloHitVector[curHitIdx].ringID == (int)ringID) {
	  curHitVector.push_back(caloHitVector[curHitIdx]);
	}
      }

      if (curHitVector.size() == 0) {

	// No hits, so just fill a header packet and no data packets
	std::vector<adc_t> curDataBlock;
	// Add the header packet to the DataBlock (leaving including a placeholder for
	// the number of packets in the DataBlock);
	adc_t null_adc = 0;
	// First 16 bits of header (reserved values)
	curDataBlock.push_back(null_adc);
	// Second 16 bits of header (ROC ID, packet type, and ring ID):
	adc_t curROCID = rocID; // 4 bit ROC ID
	adc_t headerPacketType = 5; // 4 bit Data packet header type is 5
	headerPacketType <<= 4; // Shift left by 4
	adc_t curRingID = ringID; // 3 bit ring ID
	curRingID <<= 8; // Shift left by 8
	adc_t secondEntry = (curROCID | headerPacketType | curRingID);
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
	
	// Fill in the byte count field of the header packet
	adc_t numBytes = 16; // Just the header packet
	curDataBlock[0] = numBytes;
	
	// Create mu2e::DataBlock and add to the collection
	DataBlock theBlock(DataBlock::CAL, evt.id(), dtcID, curDataBlock);
	dtcPackets->push_back(theBlock);
	
      } else {
	for (size_t curHitIdx = 0; curHitIdx < curHitVector.size(); curHitIdx++) {
	  // Generate a DataBlock for the current hit
	  
	  calhit curHit = curHitVector[curHitIdx];
	  
	  std::vector<adc_t> curDataBlock;
	  // Add the header packet to the DataBlock (leaving including a placeholder for
	  // the number of packets in the DataBlock);
	  adc_t null_adc = 0;
	  // First 16 bits of header (reserved values)
	  curDataBlock.push_back(null_adc);
	  // Second 16 bits of header (ROC ID, packet type, and ring ID):
	  adc_t curROCID = rocID; // 4 bit ROC ID
	  adc_t headerPacketType = 5; // 4 bit Data packet header type is 5
	  headerPacketType <<= 4; // Shift left by 4
	  adc_t curRingID = ringID; // 3 bit ring ID
	  curRingID <<= 8; // Shift left by 8
	  adc_t secondEntry = (curROCID | headerPacketType | curRingID);
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
	  
	  // Create a vector of adc_t values corresponding to
	  // the content of CAL data packets.
	  std::vector<adc_t> packetVector;
	  
	  // Fill the data packets:
	  // Assume the 0th apd is always read out before the second
	  adc_t crystalID = curHit.crystalId;
	  adc_t apdID = curHit.recoDigiId % 2;
	  adc_t IDNum = ((apdID << 12) | crystalID);
	  
	  packetVector.push_back(IDNum);
	  packetVector.push_back((adc_t)(curHit.recoDigiT0));

	  if(_enableDIRACEmulation==0) {
	    packetVector.push_back((adc_t)(curHit.recoDigiSamples));
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
	    packetVector.push_back( (0xFF00 & (adc_t)(maxIdx<<8)) | (0x00FF & (adc_t)(curHit.recoDigiSamples)) );
	  }


	  for (auto sampleIdx = 0; sampleIdx < curHit.recoDigiSamples; sampleIdx++) {
	      adc_t scaledVal = static_cast<adc_t>(curHit.waveform[sampleIdx]);
	      packetVector.push_back(scaledVal);
	  }

	  // Pad any empty space in the last packet with 0s
	  size_t padding_slots = 8 - ((curHit.recoDigiSamples - 5) % 8);
	  if (padding_slots < 8) {
	    for (size_t i = 0; i < padding_slots; i++) {
	      packetVector.push_back((adc_t)0);
	    }
	  }
	  
	  // Fill in the number of data packets entry in the header packet
	  adc_t numDataPackets = static_cast<adc_t>(packetVector.size() / 8);
	  curDataBlock[2] = numDataPackets;
	  
	  // Fill in the byte count field of the header packet
	  adc_t numBytes = (numDataPackets + 1) * 16;
	  curDataBlock[0] = numBytes;
	  
	  // Append the data packets after the header packet in the DataBlock
	  curDataBlock.insert(curDataBlock.end(), packetVector.begin(), packetVector.end());
	  
	  // Create mu2e::DataBlock and add to the collection
	  DataBlock theBlock(DataBlock::CAL, evt.id(), dtcID, curDataBlock);
	  dtcPackets->push_back(theBlock);
	} // Done looping over hits for this roc/ring pair

      }
	

    } // Done looping of roc/ring pairs

    evt.put(move(dtcPackets));

  } // end of ::produce

}


using mu2e::CaloPacketProducer;
DEFINE_ART_MODULE(CaloPacketProducer);
