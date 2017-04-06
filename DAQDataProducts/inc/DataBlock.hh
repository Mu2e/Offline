#ifndef DAQDataProducts_DataBlock_hh
#define DAQDataProducts_DataBlock_hh
//
// Basic storage class for a 128 bit DTC packet
//
//
// Original author Tomonari Miyashita
//

// C++ includes
#include <iostream>
#include <vector>

#include "canvas/Persistency/Provenance/EventID.h"

#include <stdexcept>

namespace mu2e {

  class DataBlock{

  public:

    enum SYSID { TRK, CAL, CRV };

    typedef uint16_t adc_t;
    typedef uint64_t timestamp;
    typedef uint8_t  dtc_id;

    // Constructors
    DataBlock() { _theID = 255; _theEventID = 0; _theSYSID=TRK; }
    
    DataBlock(const DataBlock& other) : 
      _theDataBlock(other._theDataBlock),
      _theID(other._theID),
      _theEventID(other._theEventID),
      _theSYSID(other._theSYSID)
    {}
    
    DataBlock(const DataBlock::SYSID sysid, 
	      const DataBlock::timestamp evtId, 
	      const DataBlock::dtc_id id, 
	      const std::vector<DataBlock::adc_t>& packets) {
      if(packets.size()%8==0) {
	_theDataBlock = packets;
	_theSYSID = sysid;
	_theEventID = evtId;
	_theID = id;
      } else {
	throw std::invalid_argument( "constructor received vector with incompatible length" );
      }
    }

    DataBlock(const DataBlock::SYSID sysid, 
	      const art::EventID evtId,
	      const DataBlock::dtc_id id, 
	      const std::vector<DataBlock::adc_t>& packets) {
      if(packets.size()%8==0) {
	_theDataBlock = packets;
	_theSYSID = sysid;
	_theEventID = generateUniqueID(evtId);
	_theID = id;
      } else {
	throw std::invalid_argument( "constructor received vector with incompatible length" );
      }
    }
    
    void print( std::ostream& ost = std::cout, bool doEndl = true ) const;

    size_t size() const;

    DataBlock::timestamp generateUniqueID(art::EventID evtId) const;

    std::vector<DataBlock::adc_t> timestampVector() const;

    DataBlock::timestamp getTimestamp() const;

    DataBlock::timestamp getEventID() const;

    DataBlock::dtc_id getDTCID() const;

    DataBlock::SYSID getSYSID() const;

    void setTimestamp(DataBlock::timestamp ts);

    void setDTCID(DataBlock::dtc_id id);

    void setSYSID(DataBlock::SYSID id);

    void clear();

    void appendPackets(const std::vector<DataBlock::adc_t>& _packets);

    DataBlock::adc_t at(size_t idx) const;

    DataBlock::adc_t operator[](std::size_t idx) const;

  protected:
    
    std::vector<DataBlock::adc_t> _theDataBlock;

    // ID (within the subsystem) of the DTC
    DataBlock::dtc_id _theID;

    // Unique identifier for the event (64 bits used)
    DataBlock::timestamp _theEventID;

    // Enumeration for the subsystem
    DataBlock::SYSID _theSYSID;

  };

  inline std::ostream& operator<<( std::ostream& ost, DataBlock const& packet) {
    packet.print(ost,false);
    return ost;
  }
  
} // namespace mu2e

#endif /* DAQDataProducts_DataBlock_hh */
