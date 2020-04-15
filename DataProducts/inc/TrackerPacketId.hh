#ifndef DataProducts_TrackerPacketId_hh
#define DataProducts_TrackerPacketId_hh
//
// Identifier of one straw in a tracker.
// Original author Rob Kutschke
// Re-implemented as integer bitfields by Dave Brown (LBNL)
//
#include "DataProducts/inc/StrawId.hh"
#include "cetlib_except/exception.h"

namespace mu2e {

  class TrackerPacketId{

    private:
      //  data member is a short
      uint16_t _pid;

      // define the bit field shifts and masks
    public:
      constexpr static uint16_t _channelmsk = 0x7F; // mask for panel field
      constexpr static uint16_t _channelsft = 0; // shift for preamp field
      constexpr static uint16_t _rocmsk = 0xFF00; // mask for panel field
      constexpr static uint16_t _rocsft = 8; // shift for preamp field
      constexpr static uint16_t _invalid = 0xFFFF; // invalid identifier

      constexpr static uint16_t _nchannels = StrawId::_nstraws; 
      constexpr static uint16_t _nrocs = StrawId::_nupanels;

      static bool validChannel(uint16_t ichannel) { return ichannel < _nchannels; }
      static bool validROC(uint16_t iroc) { return iroc < _nrocs; }

      TrackerPacketId(): _pid(_invalid) {}

      // construct from fields
      TrackerPacketId( uint16_t channel, uint16_t roc) : _pid(0) {
        setChannel(channel);
        setROC(roc);
      };

      // No automatic conversion of uint16_t to TrackerPacketId.
      explicit TrackerPacketId(uint16_t pid):
        _pid(pid){
        valid();
      }

      // Use compiler-generated copy c'tor, copy assignment, and d'tor.

      // test validity
      bool valid() const { return validChannel(channel()) &&
	validROC(roc()); }

      // various accessors
      uint16_t asUint16() const { return _pid;}

      uint16_t channel() const{
	return (_pid & _channelmsk) >> _channelsft;
      }

      uint16_t roc() const{
	return (_pid & _rocmsk) >> _rocsft;
      }

    private:
      // fill fields
      void setChannel(uint16_t ichannel) {
        if(validChannel(ichannel))
          _pid |= (ichannel << _channelsft) & _channelmsk;
        else
          throw cet::exception("CONFIG") << "invalid channel " << ichannel << "\n";
      }
      void setROC(uint16_t iroc) {
        if(validROC(iroc))
          _pid |= (iroc << _rocsft) & _rocmsk;
        else
          throw cet::exception("CONFIG") << "invalid roc " << iroc << "\n";
      }
  };

}
#endif /* DataProducts_StrawId_hh */
