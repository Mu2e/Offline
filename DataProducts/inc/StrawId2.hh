#ifndef TrackerGeom_StrawId2_hh
#define TrackerGeom_StrawId2_hh
//
// Identifier of one straw in a tracker.
//

//
// $Id: StrawId2.hh,v 1.7 2011/05/20 19:18:45 wb Exp $
// $Author: wb $
// $Date: 2011/05/20 19:18:45 $
//
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
      unsigned short _sid;
      // define the bit field shifts and masks
      const static unsigned short _strawmsk = 0x7F; // mask for straw field
      const static unsigned short _strawsft = 0; // shift for straw field
      const static unsigned short _panelmsk = 0x380; // mask for panel field
      const static unsigned short _panelsft = 7; // shift for panel field
      const static unsigned short _planemsk = 0xFC00; // mask for plane field
      const static unsigned short _planesft = 10; // shift for plane field
      const static unsigned short _nstraws = 96; // number of straws
      const static unsigned short _npanels = 6; // number of panels
      const static unsigned short _nplanes = 36; // number of planes
    public:
      // test values
      static bool validStraw(unsigned short istraw) { return istraw < _nstraws; }
      static bool validPanel(unsigned short ipanel) { return ipanel < _npanels; }
      static bool validPlane(unsigned short iplane) { return iplane < _nplanes; }
      // fill fields
      bool setStraw(unsigned short istraw);
      bool setPanel(unsigned short ipanel);
      bool setPlane(unsigned short iplane);

      StrawId2(std::string const& asstring);

      StrawId2(): _sid(0xFFFF) {}

      // construct from fields
      StrawId2( unsigned short plane,
	  unsigned short panel,
	  unsigned short straw);

      // Use compiler-generated copy c'tor, copy assignment, and d'tor.

      unsigned short getStrawId2() const { return _sid; }

      unsigned short getPlane() const{
	return (_sid & _planemsk) >> _planesft;
      }

      unsigned short getPanel() const{
	return (_sid & _panelmsk) >> _panelsft;
      }

      unsigned short getStraw() const{
	return (_sid & _strawmsk) >> _strawsft;
      }

      unsigned short getLayer() const{
	return _sid % 2 == 0 ? 0 : 1;
      }

      unsigned short getStation() const{
	return floor(getPlane()/2);
      }

      bool operator==( StrawId2 const& rhs) const{
	return ( _sid == rhs._sid );
      }

      bool operator!=( StrawId2 const& rhs) const{
	return !( *this == rhs);
      }

      friend std::ostream& operator<<(std::ostream& ost,
	  const StrawId2& s ){
	ost << "Plane " << s.getPlane() << " "
	  << "Panel " << s.getPanel() << " "
	  << "Straw " << s.getStraw() << std::endl;
	  return ost;
      }
  };
  inline std::ostream& operator<<(std::ostream& ost,
      const StrawId2& s );

}
#endif /* TrackerGeom_StrawId2_hh */
