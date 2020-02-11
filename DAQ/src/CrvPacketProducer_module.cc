//
// EDProducer module for converting straw hit digis into DTC formatted
// packets (which are stored in DataBlock data products)
//
// IMPORTANT NOTE: This is currently a skeleton placeholder and will need
// to be updated once the CRV and DAQ groups settle on an updated CRV packet
// specification
//
// Original author Tomonari Miyashita
//
//

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>

#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"

#include "RecoDataProducts/inc/CrvDigiCollection.hh"

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
  class CrvPacketProducer : public art::EDProducer {
  public:

    using adc_t = mu2e::DataBlock::adc_t;
    using dtc_id = mu2e::DataBlock::dtc_id;

    struct crvhit {
      DataBlock::timestamp evt;
      adc_t sipmID;
      adc_t time;
      adc_t recoDigiSamples;
      std::vector<adc_t> waveform;
      
      dtc_id dtcID;
      adc_t rocID;
    };

    explicit CrvPacketProducer(fhicl::ParameterSet const& pset);

    virtual void beginJob() override;

    virtual void endJob();

    void produce( art::Event & ) override;

  private:

    string                _outputFile;
    ofstream              outputStream;

    // For now, SiPM IDs start at 0 and end at number_of_rocs*number_of_sipms_per_roc-1
    const size_t number_of_rocs = 14;
    const size_t number_of_sipms_per_roc = 2048; // 24*64;
    const size_t number_of_rocs_per_dtc = 6;

    // 6 rocs per DTC => 3 DTCs
    // 14 rocs * 2048 sipms per roc => 28,672 sipms max

    int _generateTextFile;

    // Diagnostics level.
    int _diagLevel;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    // Label of the module that made the hits.
    std::string _makerModuleLabel;

  };

  CrvPacketProducer::CrvPacketProducer(fhicl::ParameterSet const& pset):
    art::EDProducer{pset},
    _outputFile                    (pset.get<string>("outputFile","artdaq_crv.txt")),
    _generateTextFile(pset.get<int>("generateTextFile",0)),
    _diagLevel(pset.get<int>("diagLevel",0)),
    _maxFullPrint(pset.get<int>("maxFullPrint",5)),
    _makerModuleLabel(pset.get<std::string>("makerModuleLabel","CrvDigi")) {
    produces<DataBlockCollection>();

    if(_generateTextFile>0) {
      outputStream.open(_outputFile, std::ofstream::out | std::ofstream::trunc);
    }
  }

  void CrvPacketProducer::beginJob(){

    if ( _diagLevel > 0 ) {
      cout << "CrvPacketProducer Diaglevel: "
           << _diagLevel << " "
           << _maxFullPrint
           << endl;
    }

  }

  void CrvPacketProducer::endJob(){
    if(_generateTextFile>0) {
      outputStream << flush;
      outputStream.close();
    }
  }

  void CrvPacketProducer::produce(art::Event & evt) {

    unique_ptr<DataBlockCollection> dtcPackets(new DataBlockCollection);

    bool gblresult;
    art::Handle<CrvDigiCollection> crvDigisHandle;
    gblresult = evt.getByLabel(_makerModuleLabel,"", crvDigisHandle);
    ( _diagLevel > 0 ) &&
      std::cout
      << __func__ << " getting data by getByLabel: label, instance, result " << std::endl
      << " CrvDigiCollection             _makerModuleLabel: " << _makerModuleLabel << ", "
      << ", " << gblresult << std::endl;

    if (!gblresult) throw cet::exception("DATA") << " Missing digi data";

    CrvDigiCollection const& hits_CRV = gblresult ? *crvDigisHandle : CrvDigiCollection();

    auto eventNum = evt.id().event();
    if ( _diagLevel > 2 ) {
      cout << "CrvPacketProducer: eventNum: " << eventNum
	   << " Total number of crv hit digis = " << hits_CRV.size() << endl;
    }

    std::vector<crvhit> crvHitVector; // Vector of trk hit digi data


    for(CrvDigiCollection::const_iterator iter=hits_CRV.begin(); iter!=hits_CRV.end(); iter++) {

	  const CrvDigi &crvDigi = *iter;

	  // Note: The sipmID used in the packets is a global ID incremented across all
	  // sipms while the SiPMNumber is an index from 0 to 3 for the 4 SiPMs on an
	  // individual scintillator bar
	  const CRSScintillatorBarIndex &barIndex = crvDigi.GetScintillatorBarIndex();
	  int SiPM = crvDigi.GetSiPMNumber();
	  size_t sipmID = barIndex.asInt()*4 + SiPM;

	  // 15 ROCs total
	  if(sipmID>= number_of_rocs * number_of_sipms_per_roc) {
	    throw cet::exception("DATA") << " SiPM index " << sipmID
					 << " exceeds limit of " <<  number_of_rocs
					 << " * " << number_of_sipms_per_roc
					 << " = " << number_of_rocs * number_of_sipms_per_roc;
	  }

	  // ROC ID, counting from 0, across all (for the CRV)
	  size_t globalROCID = sipmID / number_of_sipms_per_roc;
	  size_t rocID = globalROCID % number_of_rocs_per_dtc;
	  size_t dtcID = dtc_id(globalROCID / number_of_rocs_per_dtc);	  
	  
	  // Fill struct with info for current hit
	  crvhit curHit;
	  curHit.evt = eventNum;
	  curHit.sipmID = sipmID;
	  curHit.rocID = rocID;
	  curHit.dtcID = dtcID;

	  curHit.time = crvDigi.GetStartTDC();
	  curHit.recoDigiSamples = crvDigi.GetADCs().size();
	  for(size_t j = 0; j<curHit.recoDigiSamples; j++) {
	      curHit.waveform.push_back((adc_t) (crvDigi.GetADCs()[j]) );
	  }

	  crvHitVector.push_back(curHit);

	  if(_generateTextFile>0) {
	      outputStream << curHit.evt << "\t";

	      // // Added temporarily for extra debugging info
	      // outputStream << (int) curHit.dtcID << "\t";
	      // outputStream << curHit.rocID << "\t";

	      outputStream << curHit.sipmID << "\t";
	      outputStream << curHit.time << "\t";
	      //	      outputStream << curHit.recoDigiSamples << "\t";
	      for(size_t j = 0; j<curHit.waveform.size(); j++) {
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
    // // Randomize the order in which the DTCs are added to event
    // std::shuffle(dtcRocVector.begin(),dtcRocVector.end(),generator);

    // Loop over the DTC/ROC pairs and generate datablocks for each ROC
    auto targetNumROCs = dtcRocVector.size();
    for(size_t curPairNum = 0; curPairNum < targetNumROCs; curPairNum++) {

      size_t dtcID = dtcRocVector[curPairNum].first;
      size_t rocID = dtcRocVector[curPairNum].second;

      std::vector<crvhit> curHitVector;
      // Find all hits for this event coming from the specified DTC/ROC
      for (size_t curHitIdx = 0; curHitIdx < crvHitVector.size(); curHitIdx++) {
	if (crvHitVector[curHitIdx].dtcID == (unsigned int)dtcID && 
	    crvHitVector[curHitIdx].rocID == (unsigned int)rocID) {
	  curHitVector.push_back(crvHitVector[curHitIdx]);
	}
      }


      //////////////////////////////////////////////////////////////
      // Generate a DataBlock for all the hits on the current ROC
      //////////////////////////////////////////////////////////////


      // First calculate the total number of 16 bit fields and 
      // payload packets needed in the DataBlock
      
      size_t numFields = 8; // Number of 16 bit fields needed in the DataBlock
      // numFields starts at 8 because we will always have a ROC status packet
      for (size_t curHitIdx = 0; curHitIdx < curHitVector.size(); curHitIdx++) {
	size_t curNumSamples = curHitVector[curHitIdx].recoDigiSamples;
	numFields += 2 + curNumSamples/2;
	if(curNumSamples % 2 != 0) {
	  numFields += 1;
	}
	// Two additional 16 bit fields are needed for the SiPM ID, hit time, and num
	// samples for each hit
      }
      adc_t numPayloadPackets = numFields/8;
      if(numFields%8 != 0) {
	numPayloadPackets += 1;
      }
      
      
      
      std::vector<adc_t> curDataBlock;

      /////////////////////////////////////////////
      // Add the header packet to the DataBlock
      /////////////////////////////////////////////
      
      // First 16 bits of header (num bytes in DataBlock, including the header and ROC status packets)
      adc_t numBytes = (numPayloadPackets + 1) * 16;
      curDataBlock.push_back(numBytes);
      // Second 16 bits of header (ROC ID, packet type):
      adc_t curROCID = rocID; // 4 bit ROC ID
      adc_t headerPacketType = 5; // 4 bit Data packet header type is 5
      headerPacketType <<= 4; // Shift left by 4
      adc_t secondEntry = (curROCID | headerPacketType);
      secondEntry = (secondEntry | (1 << 15)); // valid bit
      curDataBlock.push_back(secondEntry);
      // Third 16 bits of header (number of data packets)
      curDataBlock.push_back(numPayloadPackets);
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
      adc_t sysID = (2 << 6) & 0x00C0;
      adc_t curDTCID = dtcID & 0x003F;
      curDataBlock.push_back(evbMode + sysID + curDTCID);
      
      
      
      /////////////////////////////////////////////////
      // Add the ROC status packet to the DataBlock
      /////////////////////////////////////////////////

      std::vector<adc_t> statusPacketVector;

      // First 16 bits contain the controller ID and packet type
      adc_t packetType = 0x0006;
      // Since the simulation doesn't currently support mapping of sipms to ROCs, we'll
      // just use the global ROC ID as the controller ID:
      adc_t controllerID = dtcID * number_of_rocs_per_dtc + rocID;
      curDataBlock.push_back(controllerID << 8 | packetType << 4);

      // Second 16 bits are the controller event word count (including the ROC status packet
      // and hits, but not the data header packet
      curDataBlock.push_back(numBytes - 16);

      // The third and fourth 16 bits contain flags indicating which FEBs are connected to
      // the ROC. This is not yet available in simulation so we'll just set them all to true:
      curDataBlock.push_back(0x00FF); // The upper 8 bits are reserved
      curDataBlock.push_back(0xFFFF);

      // The fifth 16 bits are reserved
      curDataBlock.push_back(0x0000);

      // The sixth 16 bits contain a count of how many times the CRV data has been pulled during
      // this supercycle. For now, we'll just set it to 0:
      curDataBlock.push_back(0x0000);

      // The seventh 16 bits are reserved
      curDataBlock.push_back(0x0000);

      // The eighth 16 bits contain an as yet undefined event type (we'll set it to 0) and a logical
      // or of various error bits. Since FEB errors are not implemented in the simulation, we'll set
      // all the error bits to 0:
      adc_t errors = 0x00;
      adc_t eventType = 0x00;
      curDataBlock.push_back(eventType << 8 | errors);
      

      // Append the status packets after the header packet in the DataBlock
      curDataBlock.insert(curDataBlock.end(), statusPacketVector.begin(), statusPacketVector.end());
      
      
      
      ////////////////////////////////////////////////////////////////////////////
      // Generate the rest of the DataBlock based on all the hits from this ROC
      ////////////////////////////////////////////////////////////////////////////      

      // Create a vector of adc_t values corresponding to
      // the content of CRV data payload packets.
      std::vector<adc_t> packetVector;
      for (size_t curHitIdx = 0; curHitIdx < curHitVector.size(); curHitIdx++) {
	
	crvhit curHit = curHitVector[curHitIdx];
	
	// Fill the data packets:	  
	packetVector.push_back((adc_t)(curHit.sipmID));
	
	adc_t hitTime = (adc_t)(curHit.time) & 0x03FF;
	adc_t hitSamples = (adc_t)(curHit.recoDigiSamples) << 10;
	// Make sure the number of hits fits in 6 bits
	if (curHit.recoDigiSamples>64 - 1) throw cet::exception("DATA") 
	     << " Number of samples (" << curHit.recoDigiSamples
	     << ") is too large to fit in 6 bits" << std::endl;
	
	packetVector.push_back(hitSamples | hitTime);

	for(size_t sampleIdx = 0; sampleIdx < curHit.recoDigiSamples; sampleIdx+=2) {
	  adc_t scaledVal0 = static_cast<adc_t>(curHit.waveform[sampleIdx]) & 0x00FF;
	  adc_t scaledVal1 = 0x0000;
	  if(curHit.recoDigiSamples%2==0 || sampleIdx+1 < curHit.recoDigiSamples) {
	    scaledVal1 = static_cast<adc_t>(curHit.waveform[sampleIdx+1]) << 8;
	  }
	  packetVector.push_back(scaledVal1 | scaledVal0);
	}
	
      } // Done looping over hits for this DTC/ROC pair
      
	// Pad any empty space in the last packet with 0s
      size_t padding_slots = 8 - (numFields % 8);
      if (padding_slots < 8) {
	for (size_t i = 0; i < padding_slots; i++) {
	  packetVector.push_back((adc_t)0);
	}
      }
      
      //// Fill in the number of data packets entry in the header packet
      //adc_t numDataPackets = static_cast<adc_t>(packetVector.size() / 8);
      //curDataBlock[2] = numDataPackets;
      
      //// Fill in the byte count field of the header packet
      //adc_t numBytes = (numDataPackets + 1) * 16;
      //curDataBlock[0] = numBytes;
      
      // Append the data packets after the header packet in the DataBlock
      curDataBlock.insert(curDataBlock.end(), packetVector.begin(), packetVector.end());
      
      // Create mu2e::DataBlock and add to the collection
      //	  DataBlock theBlock(DataBlock::CAL, uniqueID, dtcID, curDataBlock);
      DataBlock theBlock(DataBlock::CRV, evt.id(), dtcID, curDataBlock);
      dtcPackets->push_back(theBlock);


    } // Done looping of DTC/ROC pairs

    evt.put(move(dtcPackets));

  } // end of ::produce

}


using mu2e::CrvPacketProducer;
DEFINE_ART_MODULE(CrvPacketProducer);
