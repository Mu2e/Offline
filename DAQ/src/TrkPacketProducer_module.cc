//
// EDProducer module for converting straw hit digis into DTC formatted
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
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawDigiCollection.hh"
#include "DAQDataProducts/inc/DataBlockCollection.hh"

#include "SeedService/inc/SeedService.hh"

#include <fstream>
#include <stdexcept>

using namespace std;

namespace mu2e {

  //--------------------------------------------------------------------
  //
  //
  class TrkPacketProducer : public art::EDProducer {
  public:

    using adc_t = mu2e::DataBlock::adc_t;
    using dtc_id = mu2e::DataBlock::dtc_id;

    struct trkhit {
      DataBlock::timestamp evt;
      int strawIdx;
      unsigned long  recoDigiT0;
      unsigned long  recoDigiT1;
      unsigned long  recoDigiToT1;
      unsigned long  recoDigiToT2;
      int recoDigiSamples;
      std::vector<adc_t> waveform;
      
      int rocID;
      int ringID;

      dtc_id dtcID;
    };

    explicit TrkPacketProducer(fhicl::ParameterSet const& pset);

    virtual void beginJob() override;

    void produce( art::Event & ) override;

  private:

    string                _outputFile;
    ofstream              outputStream;

    // For now, crystal IDs start at 0 and end at number_of_rings*rocs_per_ring-1
    const size_t number_of_straws_per_roc = 96; // Each panel in the tracker has 96 straws
    const size_t numADCSamples = 12;

    // 6 optical transceivers per DTC
    // 2 optical links per ring (readout both ends)
    // => 3 Rings per DTC/server
    //
    // 12 ROCs/server
    // 3 Rings per server
    // => 4 ROCs per ring
    //
    // 20 Servers
    // 3 rings per server
    // => 60 rings
    //
    // 20 servers
    // 12 ROCs/server
    // 96 straws per ROC
    // => MAX 23040 Straws
    //
    const size_t number_of_rings = 60;
    const size_t rings_per_dtc = 3;
    const size_t rocs_per_ring = 4;

    int _generateTextFile;

    // Diagnostics level.
    int _diagLevel;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    // Label of the module that made the hits.
    std::string _makerModuleLabel;

  };

  TrkPacketProducer::TrkPacketProducer(fhicl::ParameterSet const& pset):
    _outputFile                    (pset.get<string>("outputFile","artdaq_trk.txt")),
    _generateTextFile(pset.get<int>("generateTextFile",0)),
    _diagLevel(pset.get<int>("diagLevel",0)),
    _maxFullPrint(pset.get<int>("maxFullPrint",5)),
    _makerModuleLabel(pset.get<std::string>("makerModuleLabel","makeSD")) {

    produces<DataBlockCollection>();

    if(_generateTextFile>0) {
      outputStream.open(_outputFile, std::ofstream::out | std::ofstream::trunc);
    }
  }

  void TrkPacketProducer::beginJob(){

    if ( _diagLevel > 0 ) {
      cout << "TrkPacketProducer Diaglevel: "
           << _diagLevel << " "
           << _maxFullPrint
           << endl;
    }

  }

  void TrkPacketProducer::produce(art::Event & evt) {

    unique_ptr<DataBlockCollection> dtcPackets(new DataBlockCollection);

    bool gblresult;
    art::Handle<StrawDigiCollection> strawDigisHandle;
    gblresult = evt.getByLabel(_makerModuleLabel, strawDigisHandle);
    ( _diagLevel > 0 ) && 
      std::cout 
      << __func__ << " getting data by getByLabel: label, instance, result " << std::endl
      << " StrawHitCollection             _makerModuleLabel: " << _makerModuleLabel << ", "
      << ", " << gblresult << std::endl;

    if (!gblresult) throw cet::exception("DATA") << " Missing digi data";

    StrawDigiCollection const& hits_SD = gblresult ? *strawDigisHandle : StrawDigiCollection();

    auto eventNum = evt.id().event();
    if ( _diagLevel > 2 ) {
      cout << "TrkPacketProducer: eventNum: " << eventNum
	   << " Total number of straw hit digis = " << hits_SD.size() << endl;
    }

    std::vector<trkhit> trkHitVector; // Vector of trk hit digi data
    for ( size_t i=0; i<hits_SD.size(); ++i ) {

      StrawDigi const& SD = gblresult ? hits_SD.at(i) : StrawDigi();
      TrkTypes::ADCWaveform const& theWaveform = SD.adcWaveform();

      // Fill struct with info for current hit
      trkhit curHit;
      curHit.evt = eventNum;
      curHit.strawIdx = SD.strawIndex().asInt();
      curHit.recoDigiT0 = SD.TDC(TrkTypes::cal);
      curHit.recoDigiT1 = SD.TDC(TrkTypes::hv);
      curHit.recoDigiToT1 = SD.TOT(TrkTypes::cal);
      curHit.recoDigiToT2 = SD.TOT(TrkTypes::hv);
      curHit.recoDigiSamples = theWaveform.size();
      for(size_t j = 0; j<theWaveform.size(); j++) {
	curHit.waveform.push_back(theWaveform[j]);
      }

      // 240 ROCs total
      if(curHit.strawIdx>= abs(number_of_rings * rocs_per_ring * number_of_straws_per_roc) ) {
	throw cet::exception("DATA") << " Straw index " << curHit.strawIdx
				     << " exceeds limit of " <<  number_of_rings << "*"
				     << rocs_per_ring << "*" << number_of_straws_per_roc
				     << "=" << number_of_rings * rocs_per_ring * number_of_straws_per_roc;
      }

      // Ring ID, counting from 0, across all (for the tracker)
      size_t globalRingID = int(curHit.strawIdx / (rocs_per_ring * number_of_straws_per_roc));

      curHit.ringID = globalRingID % rings_per_dtc;
      curHit.rocID = (curHit.strawIdx - (rocs_per_ring * number_of_straws_per_roc) * globalRingID) / number_of_straws_per_roc;  
      curHit.dtcID = dtc_id(globalRingID/rings_per_dtc);
      
      trkHitVector.push_back(curHit);

      if(_generateTextFile>0) {
	outputStream << curHit.evt << "\t";
	outputStream << curHit.strawIdx << "\t";
	outputStream << curHit.recoDigiT0 << "\t";
	outputStream << curHit.recoDigiT1 << "\t";
	outputStream << curHit.recoDigiToT1 << "\t";
	outputStream << curHit.recoDigiToT2 << "\t";
	outputStream << curHit.recoDigiSamples << "\t";
	for(size_t j = 0; j<curHit.waveform.size(); j++) {
	  outputStream << theWaveform[j];
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


      std::vector<trkhit> curHitVector;
      // Find all hits for this event coming from the specified Ring/ROC
      for (size_t curHitIdx = 0; curHitIdx < trkHitVector.size(); curHitIdx++) {
	//	if (trkHitVector[curHitIdx].rocID == (int)rocID && trkHitVector[curHitIdx].ringID == (int)ringID) {
	if (trkHitVector[curHitIdx].dtcID == (int)dtcID && 
	    trkHitVector[curHitIdx].rocID == (int)rocID && 
	    trkHitVector[curHitIdx].ringID == (int)ringID) {
	  curHitVector.push_back(trkHitVector[curHitIdx]);
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
	adc_t sysID = (0 << 6) & 0x00C0;
	adc_t curDTCID = dtcID & 0x003F;
	curDataBlock.push_back(evbMode + sysID + curDTCID);
	
	// Fill in the byte count field of the header packet
	adc_t numBytes = 16; // Just the header packet
	curDataBlock[0] = numBytes;
	
	// curDataBlockVector.push_back(curDataBlock);
	
	// Create mu2e::DataBlock and add to the collection
	//	DataBlock theBlock(DataBlock::TRK, uniqueID, dtcID, curDataBlock);
	DataBlock theBlock(DataBlock::TRK, evt.id(), dtcID, curDataBlock);
	dtcPackets->push_back(theBlock);
	
      } else {
	for (size_t curHitIdx = 0; curHitIdx < curHitVector.size(); curHitIdx++) {
	  // Generate a DataBlock for the current hit
	  
	  trkhit curHit = curHitVector[curHitIdx];
	  
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
	  adc_t sysID = (0 << 6) & 0x00C0;
	  adc_t curDTCID = dtcID & 0x003F;
	  curDataBlock.push_back(evbMode + sysID + curDTCID);
	  
	  // Create a vector of adc_t values corresponding to
	  // the content of TRK data packets.
	  std::vector<adc_t> packetVector;
	  
	  // Fill the data packets:
	  // Assume the 0th apd is always read out before the second
	  adc_t strawIndex = curHit.strawIdx;
	  
	  adc_t TDC0 = curHit.recoDigiT0 & 0xFFFF;
	  adc_t TDC1 = curHit.recoDigiT1 & 0xFFFF;
	  
      	  packetVector.push_back(strawIndex);
	  packetVector.push_back(TDC0);
	  packetVector.push_back(TDC1);


	  // Note: We only use 8 bits of each TOT value, and we could
	  // probably use only 4, though that wouldn't change the number
	  // of packets required per straw hit
	  uint32_t TOT0 = curHit.recoDigiToT1;
	  uint32_t TOT1 = curHit.recoDigiToT2;

	  adc_t TOT_Combined = (TOT1 << 8) | (TOT0 & 0x00FF);

      	  packetVector.push_back(TOT_Combined);	    

	  // Four 12-bit tracker ADC samples fit into every three slots (16 bits * 3)
	  // when we pack them tightly
	  for (int sampleIdx = 0; sampleIdx < curHit.recoDigiSamples; sampleIdx+=4){
	    adc_t sample0 = static_cast<adc_t>(curHit.waveform[sampleIdx]);
	    adc_t sample1 = (sampleIdx+1<curHit.recoDigiSamples) ? static_cast<adc_t>(curHit.waveform[sampleIdx+1]) : 0x0000;
	    adc_t sample2 = (sampleIdx+2<curHit.recoDigiSamples) ? static_cast<adc_t>(curHit.waveform[sampleIdx+2]) : 0x0000;
	    adc_t sample3 = (sampleIdx+3<curHit.recoDigiSamples) ? static_cast<adc_t>(curHit.waveform[sampleIdx+3]) : 0x0000;
	    
	    packetVector.push_back((sample1 << 12) | (sample0 & 0x0FFF)      );
	    packetVector.push_back((sample2 << 8) | ((sample1 >> 4) & 0x00FF));
	    packetVector.push_back((sample3 << 4) | ((sample2 >> 8) & 0x000F));
	  }

	  // Pad any empty space in the last packet with 0s
	  size_t padding_slots = 8 - (packetVector.size() % 8);
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
	  //	curDataBlockVector.push_back(curDataBlock);
	  
	  // Create mu2e::DataBlock and add to the collection
	  //	  DataBlock theBlock(DataBlock::TRK, uniqueID, dtcID, curDataBlock);
	  DataBlock theBlock(DataBlock::TRK, evt.id(), dtcID, curDataBlock);
	  dtcPackets->push_back(theBlock);
	} // Done looping over hits for this roc/ring pair

      }
	

    } // Done looping of roc/ring pairs

    evt.put(move(dtcPackets));

  } // end of ::produce

}


using mu2e::TrkPacketProducer;
DEFINE_ART_MODULE(TrkPacketProducer);
