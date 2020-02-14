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


// Definitions and typedefs needed for compatibility with HLS codeblock
#define NUM_PRESAMPLES 4
#define START_SAMPLES 5 // 0 indexed
#define NUM_SAMPLES 15
#define LOWER_TDC 16000
#define UPPER_TDC 64000

typedef uint16_t tdc_type;
typedef uint8_t  tot_type;
typedef uint16_t adc_type;
typedef uint16_t calib_constant_type;
typedef uint8_t  flag_mask_type;


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
      int strawID;
      unsigned long  recoDigiT0;
      unsigned long  recoDigiT1;
      unsigned long  recoDigiToT1;
      unsigned long  recoDigiToT2;
      size_t recoDigiSamples;
      std::vector<adc_t> waveform;

      dtc_id dtcID;      
      int rocID;      

      flag_mask_type flags;
    };

    explicit TrkPacketProducer(fhicl::ParameterSet const& pset);

    virtual void beginJob() override;

    virtual void endJob();

    void produce( art::Event & ) override;

    flag_mask_type filter( tdc_type tdc0, tdc_type tdc1,
			   tot_type tot0, tot_type tot1,
			   adc_type adc[NUM_SAMPLES],
			   
			   calib_constant_type clockstart, 
			   calib_constant_type panelTDCoffset, calib_constant_type hvoffset, calib_constant_type caloffset,
			   calib_constant_type energy_max_LSHIFT8, calib_constant_type energy_min_LSHIFT8,
			   calib_constant_type gain_RSHIFT15,
			   calib_constant_type inverse_ionization_energy_LSHIFT26);

  private:

    string                _outputFile;
    ofstream              outputStream;

    //    const size_t number_of_rocs = 216;
    const size_t number_of_rocs = 240;
    const size_t number_of_straws_per_roc = 96; // Each panel in the tracker has 96 straws
    const size_t number_of_rocs_per_dtc = 6;
    const size_t numADCSamples = 15;

    // 96 straws per panel
    // 1 ROC per panel
    // 216 panels
    //
    // 6 ROCs per DTC
    // 36 DTCs

    int _generateTextFile;

    // Diagnostics level.
    int _diagLevel;

    int _enableFPGAEmulation;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    // Label of the module that made the hits.
    std::string _makerModuleLabel;

  };

  TrkPacketProducer::TrkPacketProducer(fhicl::ParameterSet const& pset):
    art::EDProducer{pset},
    _outputFile                    (pset.get<string>("outputFile","artdaq_trk.txt")),
    _generateTextFile(pset.get<int>("generateTextFile",0)),
    _diagLevel(pset.get<int>("diagLevel",0)),
    _enableFPGAEmulation(pset.get<int>("enableFPGAEmulation",0)),
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

      if(_enableFPGAEmulation) {
	cout << "FPGA Emulation Enabled" << endl;
      }

    }

  }

  void TrkPacketProducer::endJob(){
    if(_generateTextFile>0) {
      outputStream << flush;
      outputStream.close();
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
      curHit.strawID = SD.strawId().asUint16();      
      curHit.recoDigiT0 = SD.TDC(StrawEnd::cal);
      curHit.recoDigiT1 = SD.TDC(StrawEnd::hv);
      curHit.recoDigiToT1 = SD.TOT(StrawEnd::cal);
      curHit.recoDigiToT2 = SD.TOT(StrawEnd::hv);
      //      curHit.recoDigiSamples = theWaveform.size();
      // Since the TRK digis still contain 16 samples, this
      // must be hardcoded to 15 until the digi code is updated.
      curHit.recoDigiSamples = numADCSamples;
      for(size_t j = 0; j<curHit.recoDigiSamples; j++) {
	curHit.waveform.push_back(theWaveform[j]);
      }

      // 96 straws per ROC/panel
      // 6 panels / plane
      //      // 36 planes => 216 panels/ROCs
      // There are actually 40 planes in the input file => 240 panels/ROCs

      int panel = SD.strawId().getPanel();
      int plane = SD.strawId().getPlane();

      // ROC ID, counting from 0 across all DTCs (for the tracker)
      size_t globalROCID = (plane*6) + panel;
      size_t localROCID = panel;

      curHit.rocID = localROCID;
      curHit.dtcID = dtc_id(globalROCID/number_of_rocs_per_dtc);

//      std::cout << "GREPME GLOBALROCID "  << globalROCID
//		<<  " ROCID "             << curHit.rocID  << "/" << number_of_rocs_per_dtc
//		<<  " DTCID "             << (int)curHit.dtcID  << "/" << (int)(number_of_rocs/number_of_rocs_per_dtc-1)
//		<<  " panelID "           << panel
//		<<  " planeID "           << plane << std::endl;

      // 240 ROCs total
      if(globalROCID >= number_of_rocs) {
	throw cet::exception("DATA") << " Global ROC ID " << globalROCID
				     << " exceeds limit of " << number_of_rocs;
      }
       
      trkHitVector.push_back(curHit);

//      if(_generateTextFile>0) {
//	outputStream << curHit.evt << "\t";
//	outputStream << curHit.strawID << "\t";
//	outputStream << curHit.recoDigiT0 << "\t";
//	outputStream << curHit.recoDigiT1 << "\t";
//	outputStream << curHit.recoDigiToT1 << "\t";
//	outputStream << curHit.recoDigiToT2 << "\t";
//	outputStream << curHit.recoDigiSamples << "\t";
//	for(size_t j = 0; j<curHit.waveform.size(); j++) {
//	  outputStream << theWaveform[j];
//	  if(j<curHit.waveform.size()-1) {
//	    outputStream << "\t";
//	  }
//	}
//
//	if(_enableFPGAEmulation) {
//	  outputStream << unsigned(curHit.flags);
//	}
//	
//	outputStream << endl;
//      }



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

      std::vector<trkhit> curHitVector;
      // Find all hits for this event coming from the specified DTC/ROC combination
      for (size_t curHitIdx = 0; curHitIdx < trkHitVector.size(); curHitIdx++) {
	if (trkHitVector[curHitIdx].dtcID == (int)dtcID && 
	    trkHitVector[curHitIdx].rocID == (int)rocID ) {
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
	  adc_t sysID = (0 << 6) & 0x00C0;
	  adc_t curDTCID = dtcID & 0x003F;
	  curDataBlock.push_back(evbMode + sysID + curDTCID);
	  
	  // Create a vector of adc_t values corresponding to
	  // the content of TRK data packets.
	  std::vector<adc_t> packetVector;
	  
	  // Fill the data packets:
	  // Assume the 0th apd is always read out before the second
	  uint16_t strawID = curHit.strawID;
	  
	  adc_t TDC0 = curHit.recoDigiT0 & 0xFFFF;
	  adc_t TDC1 = curHit.recoDigiT1 & 0xFFFF;
	  
      	  packetVector.push_back(strawID);
	  packetVector.push_back(TDC0);
	  packetVector.push_back(TDC1);


	  // Note: We only use 8 bits of each TOT value, and we could
	  // probably use only 4, though that wouldn't change the number
	  // of packets required per straw hit
	  uint8_t TOT0 = curHit.recoDigiToT1;
	  uint8_t TOT1 = curHit.recoDigiToT2;

	  adc_t TOT_Combined = (TOT1 << 8) | (TOT0 & 0x00FF);

      	  packetVector.push_back(TOT_Combined);	    

	  // Four 12-bit tracker ADC samples fit into every three slots (16 bits * 3)
	  // when we pack them tightly
	  for (size_t sampleIdx = 0; sampleIdx < curHit.recoDigiSamples; sampleIdx+=4){
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
	  
	  if(_enableFPGAEmulation) {
	    
	    adc_type adc[NUM_SAMPLES];
	    for(size_t i=0; i<NUM_SAMPLES; i++) {
	      adc[i] = curHit.waveform[i];
	    }

	    // Note: Eventually, there will be a calibration database to provide these sorts of values
	    // For now, these are just placeholders
	    calib_constant_type clockstart = 320; //10 / .03125 (tdclsb)
	    calib_constant_type maxenergyls8 = 583, minenergyls8 = 0;
	    calib_constant_type gainrs15 = 1389;
	    calib_constant_type inverseionizationenergyls26 = 633;
	    calib_constant_type panelTDCoffset = 0,  hvoffset = 0,  caloffset = 0;

	    flag_mask_type f = filter(TDC0,TDC1,TOT0,TOT1,adc,clockstart,panelTDCoffset,hvoffset,caloffset,maxenergyls8,minenergyls8,gainrs15,inverseionizationenergyls26);
	    curHit.flags = f;

	    adc_t preprocessing_flags = 0x0000 | (f << 8);

	    packetVector[packetVector.size()-1] = preprocessing_flags | packetVector[packetVector.size()-1];
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


	  if(_generateTextFile>0) {
	    outputStream << curHit.evt << "\t";
	    outputStream << curHit.strawID << "\t";
	    outputStream << curHit.recoDigiT0 << "\t";
	    outputStream << curHit.recoDigiT1 << "\t";
	    outputStream << curHit.recoDigiToT1 << "\t";
	    outputStream << curHit.recoDigiToT2 << "\t";
	    outputStream << curHit.recoDigiSamples << "\t";
	    for(size_t j = 0; j<curHit.waveform.size(); j++) {
	      outputStream << curHit.waveform[j];
	      if(j<curHit.waveform.size()-1) {
		outputStream << "\t";
	      }
	    }
	    if(_enableFPGAEmulation) {
	      outputStream << "\t" << unsigned(curHit.flags);
	    }
	    outputStream << endl;
	  }

	} // Done looping over hits for this DTC/ROC pair

      }
	

    } // Done looping of DTC/ROC pairs

    evt.put(move(dtcPackets));

  } // end of ::produce



  flag_mask_type TrkPacketProducer::filter( tdc_type tdc0, tdc_type tdc1,
					    tot_type tot0, tot_type tot1,
					    adc_type adc[NUM_SAMPLES],
					    
					    calib_constant_type clockstart, 
					    calib_constant_type panelTDCoffset, calib_constant_type hvoffset, calib_constant_type caloffset,
					    calib_constant_type energy_max_LSHIFT8, calib_constant_type energy_min_LSHIFT8,
					    calib_constant_type gain_RSHIFT15,
					    calib_constant_type inverse_ionization_energy_LSHIFT26) {

    int failed_time = 0;
    int failed_energy = 0;
    int thistdc = tdc0 < tdc1 ? tdc0 + hvoffset : tdc1 + caloffset;
    thistdc = thistdc + clockstart + panelTDCoffset;
    if (thistdc < LOWER_TDC || thistdc > UPPER_TDC) {
      failed_time = 1;
    }
    
    //sum up presamples
    int pedsum = 0;
    for (int i = 0; i < NUM_PRESAMPLES; i++) {
      pedsum += adc[i];
    }
    int pedestal = pedsum / NUM_PRESAMPLES;
    
    int peak = 0;
    for (int i = START_SAMPLES; i < NUM_SAMPLES; i++) {
      if (adc[i] > peak) {
	peak = adc[i];
      } else {
	break;
      }
    }
    
    int energy = peak - pedestal;
    
    int energy_max_adjusted = ((((energy_max_LSHIFT8 * gain_RSHIFT15) >> 9) * inverse_ionization_energy_LSHIFT26) >> 10);
    int energy_min_adjusted = ((((energy_min_LSHIFT8 * gain_RSHIFT15) >> 9) * inverse_ionization_energy_LSHIFT26) >> 10);
    if (energy > energy_max_adjusted || energy < energy_min_adjusted) {
      failed_energy = 1;
    }
    
    return (failed_energy<<1) | failed_time;
    
  }

}


using mu2e::TrkPacketProducer;
DEFINE_ART_MODULE(TrkPacketProducer);
