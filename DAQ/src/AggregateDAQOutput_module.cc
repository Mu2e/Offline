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

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
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
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "DAQDataProducts/inc/DataBlockCollection.hh"

#include "SeedService/inc/SeedService.hh"

#include <fstream>
#include <stdexcept>

using namespace std;

namespace mu2e {

  //--------------------------------------------------------------------
  //
  //
  class AggregateDAQOutput : public art::EDAnalyzer {
  public:

    explicit AggregateDAQOutput(fhicl::ParameterSet const& pset);

    virtual void beginJob() override;

    virtual void endJob();

    void analyze( art::Event const& evt) override;

  protected:

    void flushBuffer();

    std::vector<DataBlock::adc_t> generateDMABlockHeader(size_t theCount) const;

    string                _outputFile;
    ofstream              outputStream;

    size_t _maxDMABlockSize;
    // Within each event (corresponding to a unique timestamp) the DataBlocks
    // are divided into Direct Memory Access (DMA) blocks with a max size
    // in bytes corresponding to _dmaBlockSize.
    // NOTE: THE DMA BLOCK SIZE INCLUDES THE DMA BLOCK HEADER !!!

    size_t _bufferSize;
    std::vector<DataBlock::adc_t> outputBuffer;

    size_t _generateTimestampTable;

    // Table used for mapping between DTC timestamp and art EventID
    string _tableFile;
    std::vector< std::pair<DataBlock::timestamp,DataBlock::timestamp> > tsTable;


    size_t _timestampOffset;

    size_t numWordsWritten;
    size_t numEventsProcessed;

    string _TrackerDataBlockMakerModule;
    string _CalorimeterDataBlockMakerModule;
    string _CosmicRayVetoDataBlockMakerModule;

    int _includeTracker;
    int _includeCalorimeter;
    int _includeCosmicRayVeto;

    // Diagnostics level.
    int _diagLevel;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

  };

  AggregateDAQOutput::AggregateDAQOutput(fhicl::ParameterSet const& pset):
    art::EDAnalyzer(pset),
    _outputFile                     (pset.get<string>("outputFile","DTC_packets.bin")),
    _maxDMABlockSize                (pset.get<size_t>("maxDMABlockSize",32000)), // Maximum size in bytes of a DMA block
    _bufferSize                     (pset.get<size_t>("bufferSize",1000)),
    _generateTimestampTable         (pset.get<size_t>("generateTimestampTable",0)),
    _tableFile                      (pset.get<string>("tableFile","tsTable.bin")),
    _timestampOffset                (pset.get<size_t>("timestampOffset",0)),
    numWordsWritten(0),
    numEventsProcessed(0),
    _TrackerDataBlockMakerModule    (pset.get<string>("TrackerDataBlockMakerModule","TrackerPacketProducer")),
    _CalorimeterDataBlockMakerModule(pset.get<string>("CalorimeterDataBlockMakerModule","CalorimeterPacketProducer")),
    _CosmicRayVetoDataBlockMakerModule(pset.get<string>("CosmicRayVetoDataBlockMakerModule","CosmicRayVetoPacketProducer")),
    _includeTracker(pset.get<int>("includeTracker",1)),
    _includeCalorimeter(pset.get<int>("includeCalorimeter",1)),
    _includeCosmicRayVeto(pset.get<int>("includeCosmicRayVeto",0)),
    _diagLevel(pset.get<int>("diagLevel",0)),
    _maxFullPrint(pset.get<int>("maxFullPrint",5))
  {
    outputStream.open(_outputFile, std::ios::out | std::ios::binary);
  }

  void AggregateDAQOutput::beginJob(){
    if ( _diagLevel > 0 ) {
      cout << "AggregateDAQOutput Diaglevel: "
           << _diagLevel << " "
           << _maxFullPrint
           << endl;
    }
  }

  void AggregateDAQOutput::endJob(){
    flushBuffer();
    outputStream.close();

    if(_generateTimestampTable) {
      ofstream tsTableStream;
      tsTableStream.open(_tableFile, std::ios::out | std::ios::binary);
      for(size_t idx=0; idx<tsTable.size(); idx++) {
	tsTableStream.write(reinterpret_cast<const char *>(&(tsTable[idx].first)), sizeof(DataBlock::timestamp));
	tsTableStream.write(reinterpret_cast<const char *>(&(tsTable[idx].second)), sizeof(DataBlock::timestamp));

	if (_diagLevel > 3) {
	  cout << "TIMESTAMP_MAPPING: timestamp: "
	       << tsTable[idx].first
	       << " uniqueid: "
	       << tsTable[idx].second
	       << endl;
	}

      }
      tsTableStream << std::flush;
      tsTableStream.close();
    }

    if (_diagLevel > 0) {
      cout << "AggregateDAQOutput: "
	   << "Finished writing "
	   << numWordsWritten
	   << " words from "
	   << numEventsProcessed
	   << " events to "
	   << _outputFile
	   << endl;
    }

  }

  std::vector<DataBlock::adc_t> AggregateDAQOutput::generateDMABlockHeader(size_t theCount) const {

    uint64_t byteCount = theCount;

    std::vector<DataBlock::adc_t> header;
    header.push_back(static_cast<DataBlock::adc_t>( byteCount        & 0xFFFF));
    header.push_back(static_cast<DataBlock::adc_t>((byteCount >> 16) & 0xFFFF));
    header.push_back(static_cast<DataBlock::adc_t>((byteCount >> 32) & 0xFFFF));
    header.push_back(static_cast<DataBlock::adc_t>((byteCount >> 48) & 0xFFFF));
    
    return header;
  }

  void AggregateDAQOutput::flushBuffer() {

    for(size_t idx = 0; idx<outputBuffer.size(); idx++) {
      outputStream.write(reinterpret_cast<const char *>(&(outputBuffer[idx])), sizeof(DataBlock::adc_t));
    }
    outputStream << std::flush;

    numWordsWritten+=outputBuffer.size()*2;
    outputBuffer.clear();
  }

  void AggregateDAQOutput::analyze(art::Event const& evt) {
    
    bool tsTableEntryRecorded = false;

    DataBlock::timestamp ts = numEventsProcessed + _timestampOffset;

    std::vector<DataBlockCollection> theCollections;

    if(_includeTracker>0) {
      bool gblTresult;
      art::Handle<DataBlockCollection> trkDataBlockHandle;
      gblTresult = evt.getByLabel(_TrackerDataBlockMakerModule, trkDataBlockHandle);
      if (!gblTresult) throw cet::exception("DATA") << " Missing tracker DataBlock data ";
      DataBlockCollection const& trk_datablocks = gblTresult ? *trkDataBlockHandle : DataBlockCollection();

      if ( _diagLevel > 1 ) {
	cout << "AggregateDAQOutput: Total number of tracker DataBlocks = " << trk_datablocks.size() << endl;
      }
      
      theCollections.push_back(trk_datablocks);
    }

    if(_includeCalorimeter>0) {
      bool gblCresult;
      art::Handle<DataBlockCollection> caloDataBlockHandle;
      gblCresult = evt.getByLabel(_CalorimeterDataBlockMakerModule, caloDataBlockHandle);
      if (!gblCresult) throw cet::exception("DATA") << " Missing calorimeter DataBlock data ";
      DataBlockCollection const& calo_datablocks = gblCresult ? *caloDataBlockHandle : DataBlockCollection();

      if ( _diagLevel > 1 ) {
	cout << "AggregateDAQOutput: Total number of calorimeter DataBlocks = " << calo_datablocks.size() << endl;
      }

      theCollections.push_back(calo_datablocks);
    }

    if(_includeCosmicRayVeto>0) {
      bool gblVresult;
      art::Handle<DataBlockCollection> crvDataBlockHandle;
      gblVresult = evt.getByLabel(_CosmicRayVetoDataBlockMakerModule, crvDataBlockHandle);
      if (!gblVresult) throw cet::exception("DATA") << " Missing cosmic ray veto DataBlock data ";
      DataBlockCollection const& crv_datablocks = gblVresult ? *crvDataBlockHandle : DataBlockCollection();

      if ( _diagLevel > 1 ) {
	cout << "AggregateDAQOutput: Total number of cosmic ray veto DataBlocks = " << crv_datablocks.size() << endl;
      }

      theCollections.push_back(crv_datablocks);
    }

    std::vector<DataBlockCollection> collectionsByDTC;

    // Loop over the different subsystems
    for(size_t i=0; i<theCollections.size(); i++) {

      // Create a separate DataBlockCollection for each DTC and store
      // it in a dictionary where the DTC ID is the key
      std::map<DataBlock::dtc_id,DataBlockCollection> curDictionary;
      for(size_t j=0; j<theCollections[i].size(); j++) {
	curDictionary[ theCollections[i][j].getDTCID() ].push_back( theCollections[i][j] );
      }
      
      // Find all the keys in the dictionary
      std::vector<DataBlock::dtc_id> keys;
      
      for(auto it = curDictionary.begin(); it != curDictionary.end(); it++) {
	keys.push_back(it->first);
      }

      //
      // At this point, one could sort the keys in the vector if needed
      //

      for(size_t keyIdx=0; keyIdx<keys.size(); keyIdx++) {
	collectionsByDTC.push_back(curDictionary[keys[keyIdx]]);
      }

    }

    // Generate the timestamp conversion table
    for(size_t collectionIdx = 0; collectionIdx<collectionsByDTC.size(); collectionIdx++) {
      DataBlockCollection datablocks = collectionsByDTC[collectionIdx];
      if(datablocks.size() > 0 && _generateTimestampTable>0 && !tsTableEntryRecorded) {
	  std::pair<DataBlock::timestamp,DataBlock::timestamp> curPair(ts,datablocks[0].getEventID());
	  tsTable.push_back(curPair);
	  tsTableEntryRecorded = true;
      }
    }

    // Determine how to divide DataBlocks between DMABlocks within each timestamp
    std::vector<size_t> dataBlockPartition;
    std::vector<size_t> dataBlockPartitionSizes;

    size_t curDMABlockSize = 0;
    size_t numDataBlocksInCurDMABlock = 0;
    for(size_t collectionIdx = 0; collectionIdx<collectionsByDTC.size(); collectionIdx++) {
      
      DataBlockCollection datablocks = collectionsByDTC[collectionIdx];

      for(size_t dataBlockIdx = 0; dataBlockIdx<datablocks.size(); dataBlockIdx++) {
	DataBlock curDataBlock = datablocks[dataBlockIdx];

	if(numDataBlocksInCurDMABlock == 0) {
	  // Starting a new DMA Block, so allocate
	  // space for a new DMA block header
	  curDMABlockSize = 8;
	}

	numDataBlocksInCurDMABlock++; // Increment number of DataBlocks in the current DMA block
	curDMABlockSize += curDataBlock.size() * 2; // Size of current data block in 8bit words

	if(curDMABlockSize > _maxDMABlockSize) {
	  throw cet::exception("DATA") << "Current DMA Block size (" 
				       << curDMABlockSize << ") exceeds max DMA Block size ("
				       << _maxDMABlockSize << ")" << std::endl;
	}

	bool atEndOfDMABlock = false;

	if(dataBlockIdx == datablocks.size() - 1 &&
	   collectionIdx == collectionsByDTC.size() - 1) {
	  // There are no more DataBlocks, so this marks the end of
	  // the current DMA block
	  atEndOfDMABlock = true;
	} 

	if( dataBlockIdx == datablocks.size() - 1 &&
	    collectionIdx < collectionsByDTC.size() - 1 &&
	    collectionsByDTC[collectionIdx+1].size()>0 &&
	    curDMABlockSize + (2 * collectionsByDTC[collectionIdx+1][0].size()) > _maxDMABlockSize) {
	  // Adding the DataBlock would go over the size limit
	  // so this marks the end of the current DMA block
	  atEndOfDMABlock = true;
	} 	 

	if(dataBlockIdx < datablocks.size() - 1 &&
	   curDMABlockSize + (2 * datablocks[dataBlockIdx + 1].size()) > _maxDMABlockSize) {
	  // Adding the next DataBlock would put us over the limit so this is
	  // the end of the current DMA block
	  atEndOfDMABlock = true;
	}

	if(atEndOfDMABlock) {
	  dataBlockPartition.push_back(numDataBlocksInCurDMABlock);
	  dataBlockPartitionSizes.push_back(curDMABlockSize);
	  numDataBlocksInCurDMABlock = 0;
	}

      } // end loop over DataBlocks
    } // end loop over DTCs    


    // Break the DataBlocks into DMABlocks and add DMABlock headers

    // Index of the current DMA block in the partition array
    size_t curDMABlockIdx = 0;

    // Number of DataBlocks added to the current DMA block
    size_t curDataBlockCount = 0; 
    // Number of DataBlocks added to the current DMA block
    size_t curDMABlockByteCount = 0; 

    bool atBeginningOfDMABlock = true;

    std::vector<DataBlock::adc_t> masterVector;
    std::vector<DataBlock::adc_t> DMABlockVector;

    for(size_t collectionIdx = 0; collectionIdx<collectionsByDTC.size(); collectionIdx++) {
      DataBlockCollection datablocks = collectionsByDTC[collectionIdx];
      for(size_t dataBlockIdx = 0; dataBlockIdx<datablocks.size(); dataBlockIdx++) {
	DataBlock curDataBlock = datablocks[dataBlockIdx];

	if(atBeginningOfDMABlock) {
	  atBeginningOfDMABlock = false;

	  DMABlockVector.clear();
	  
	  curDataBlockCount = 0;
	    
	  std::vector<DataBlock::adc_t> header = generateDMABlockHeader( dataBlockPartitionSizes[curDMABlockIdx] );
	  DMABlockVector.insert(DMABlockVector.end(), header.begin(), header.end());

	  curDMABlockByteCount = 8;
	}

	// Add the current DataBlock to the current SuperBlock
	curDataBlock.setTimestamp(ts); // Overwrite the timestamp
	for(size_t adcNum = 0; adcNum < curDataBlock.size(); adcNum++) {
	  DMABlockVector.push_back(curDataBlock[adcNum]);
	}
	curDMABlockByteCount += curDataBlock.size() * 2; 
	curDataBlockCount++;

	if ( _diagLevel > 1 && dataBlockIdx==0) {
	  DataBlock::SYSID curSYSID = curDataBlock.getSYSID();
	  DataBlock::dtc_id curDTCID = curDataBlock.getDTCID();

	  std::cout << "================================================" << std::endl;
	  std::cout << "\t\tTimestamp: " << ts << std::endl;
	  std::cout << "\t\tDTCID: " << (int)curDTCID << std::endl;
	  std::cout << "\t\tSYSID: " << (int)curSYSID << std::endl;
	}

	if(curDataBlockCount==dataBlockPartition[curDMABlockIdx]) {
	  // Reached end of current DMA block

	  if ( _diagLevel > 1 ) {
	    std::cout << "Number of bytes in DMABlock: " << curDMABlockByteCount << std::endl;
	  }

	  masterVector.insert(masterVector.end(), DMABlockVector.begin(), DMABlockVector.end());
	  curDMABlockIdx++;

	  atBeginningOfDMABlock = true;
	}
	

      } // End loop over DataBlocks

    } // End loop over DTC collections
      


    // Write all values, including superblock header and DMA header values, to output buffer
    for ( size_t idx=0; idx<masterVector.size(); idx++ ) {
      if(outputBuffer.size()>= _bufferSize) {
	flushBuffer();
      }
      outputBuffer.push_back(masterVector[idx]);
    }

    numEventsProcessed += 1;

  } // end of ::analyze.

}


using mu2e::AggregateDAQOutput;
DEFINE_ART_MODULE(AggregateDAQOutput);
