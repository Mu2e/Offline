//
// Basic storage class for a 128 bit DTC packet
//
//
// Original author Tomonari Miyashita
//

// C++ includes
#include <ostream>

// Framework includes.
#include "cetlib_except/exception.h"

// Mu2e includes
#include "DAQDataProducts/inc/DataBlock.hh"


namespace mu2e {

  // Print the raw content of the packet.
  void DataBlock::print( std::ostream& ost, bool doEndl ) const {

    ost << "DataBlock:";
    ost << "DTCID: " << _theID;
    for(size_t i=0; i<this->size(); i++) {
      ost << " Data" << i << ": " << this->at(i);
    };
    if ( doEndl ){
      ost << std::endl;
    }

  }

  size_t DataBlock::size() const {
    return _theDataBlock.size();
  }

  DataBlock::timestamp DataBlock::generateUniqueID(art::EventID evtId) const {
    // Define unique 64-bit identifier for this event. Use 32 bits from
    // RunNumber, 14 bits from SubRunNumber, and 18 bits from EventNumber
    DataBlock::timestamp curRunID = evtId.run();
    DataBlock::timestamp curSubRunID = evtId.subRun();
    DataBlock::timestamp curEvtNum = evtId.event();
    if(curRunID > 4294967295) { throw cet::exception("DATA") << "RunNumber requires more than 32 bits: " << curRunID << " "; }
    if(curEvtNum > 262143) { throw cet::exception("DATA") << "EventNumber requires more than 18 bits: " << curEvtNum << " "; }
    if(curSubRunID > 16383) { throw cet::exception("DATA") << "SubRunNumber requires more than 14 bits: " << curSubRunID << " "; }
    DataBlock::timestamp uniqueID = curEvtNum & 0x3FFFF; // Use 18 bits for EventNum
    uniqueID = uniqueID | ( (curSubRunID & 0x3FFF ) << 18 ); // Use 14 bits for SubRun and shift to left of event num
    uniqueID = uniqueID | ( (curRunID    & 0xFFFFFFFF ) << 32 ); // Use remaining 32 bits for Run number
    return uniqueID;
  }

  std::vector<DataBlock::adc_t> DataBlock::timestampVector() const {
    std::vector<DataBlock::adc_t> theVector;
    if(this->size()<8) {
      return theVector;
    } else {
      for(size_t i=3; i<6; i++) {
	theVector.push_back(this->at(i));
      }
      return theVector;
    }
  }
  
  DataBlock::timestamp DataBlock::getTimestamp() const {
    DataBlock::timestamp ts = 0;
    std::vector<DataBlock::adc_t> theVector = this->timestampVector();
    if(theVector.size()==0) {
      return ts;
    } else {
      DataBlock::timestamp ts0 = theVector[0];
      DataBlock::timestamp ts1 = theVector[1];
      DataBlock::timestamp ts2 = theVector[2];
      ts1 <<= 16;
      ts2 <<= 32;
      ts = ts0 + ts1 + ts2;
      return ts;
    }
  }

  DataBlock::timestamp DataBlock::getEventID() const {
    return _theEventID;
  }

  DataBlock::dtc_id DataBlock::getDTCID() const {
    return _theID;
  }

  DataBlock::SYSID DataBlock::getSYSID() const {
    return _theSYSID;
  }
  
  void DataBlock::setTimestamp(DataBlock::timestamp ts) {
    _theDataBlock[3] = static_cast<adc_t>( ts        & 0xFFFF);
    _theDataBlock[4] = static_cast<adc_t>((ts >> 16) & 0xFFFF);
    _theDataBlock[5] = static_cast<adc_t>((ts >> 32) & 0xFFFF);      
  }

  void DataBlock::setDTCID(DataBlock::dtc_id id) {
    _theID = id;
  }

  void DataBlock::setSYSID(DataBlock::SYSID id) {
    _theSYSID = id;
  }

  void DataBlock::clear() {
    _theDataBlock.clear();
  }
  
  void DataBlock::appendPackets(const std::vector<DataBlock::adc_t>& _packets) {
    if(_packets.size()%8==0) {
      _theDataBlock.insert(_theDataBlock.end(), _packets.begin(), _packets.end());
    } else {
      throw std::invalid_argument( "appendPackets received vector with incompatible length" );
    }      
  }
  
  DataBlock::adc_t DataBlock::at(size_t idx) const {
    return _theDataBlock[idx];
  }
  
  DataBlock::adc_t DataBlock::operator[](std::size_t idx) const { 
    return this->at(idx); 
  }
  

} // namespace mu2e
