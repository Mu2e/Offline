#ifndef TrackerGeom_StrawId2_hh
#define TrackerGeom_StrawId2_hh
//
// Identifier of one straw in a tracker.
// Original author Rob Kutschke
// Re-implemented as integer bitfields by Dave Brown (LBNL)
//
#include <ostream>
#include <string>
#include <math.h>

namespace mu2e {

  class StrawId2{
    private:
      //  data member is a short
      uint16_t _sid;
      // define the bit field shifts and masks
      constexpr static uint16_t _strawmsk = 0x7F; // mask for straw field
      constexpr static uint16_t _panelmsk = 0x380; // mask for panel field
      constexpr static uint16_t _panelsft = 7; // shift for panel field
      constexpr static uint16_t _planemsk = 0xFC00; // mask for plane field
      constexpr static uint16_t _planesft = 10; // shift for plane field
      constexpr static uint16_t _nstraws = 96; // number of straws
      constexpr static uint16_t _npanels = 6; // number of panels
      constexpr static uint16_t _nplanes = 36; // number of planes
      constexpr static uint16_t _invalid = 0xFFFF; // invalid identifier
    
    public:
      // test values
      static bool validStraw(uint16_t istraw) { return istraw < _nstraws; }
      static bool validPanel(uint16_t ipanel) { return ipanel < _npanels; }
      static bool validPlane(uint16_t iplane) { return iplane < _nplanes; }

      StrawId2(std::string const& asstring);

      StrawId2(): _sid(_invalid) {}

      // construct from fields
      StrawId2( uint16_t plane,
	  uint16_t panel,
	  uint16_t straw);

      // Use compiler-generated copy c'tor, copy assignment, and d'tor.

    // test validity
      bool valid() const { return validStraw(straw()) &&
	validPanel(panel()) && validPlane(plane()); }

      uint16_t strawId2() const { return _sid; }

      uint16_t plane() const{
	return (_sid & _planemsk) >> _planesft;
      }

      uint16_t panel() const{
	return (_sid & _panelmsk) >> _panelsft;
      }

      uint16_t straw() const{
	return (_sid & _strawmsk);
      }

      uint16_t layer() const{
	return _sid % 2 == 0 ? 0 : 1;
      }

      uint16_t stationn() const{
	return floor(plane()/2);
      }

      bool operator==( StrawId2 const& rhs) const{
	return ( _sid == rhs._sid );
      }

      bool operator!=( StrawId2 const& rhs) const{
	return !( *this == rhs);
      }

      friend std::ostream& operator<<(std::ostream& ost,
	  const StrawId2& s ){
	ost << "StrawId with Plane " << s.plane() << " "
	  << "Panel " << s.panel() << " "
	  << "Straw " << s.straw() << std::endl;
	  return ost;
      }

    private:
      // fill fields
      void setStraw(uint16_t istraw);
      void setPanel(uint16_t ipanel);
      void setPlane(uint16_t iplane);

  };
  inline std::ostream& operator<<(std::ostream& ost,
      const StrawId2& s );

}
#endif /* TrackerGeom_StrawId2_hh */
