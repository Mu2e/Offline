#ifndef TrackerGeom_StrawId_hh
#define TrackerGeom_StrawId_hh
//
// Identifier of one straw in a tracker.
// Original author Rob Kutschke
// Re-implemented as integer bitfields by Dave Brown (LBNL)
//
#include <ostream>
#include <string>
#include <math.h>

namespace mu2e {

  class StrawId{

    friend class Tracker;
    friend class TTracker;
    friend class TTrackerMaker;
    friend class Plane;
    friend class Panel;
    friend class StrawSD;

    private:
      //  data member is a short
      uint16_t _sid;

      // define the bit field shifts and masks
      constexpr static uint16_t _strawmsk = 0x7F; // mask for straw field
      constexpr static uint16_t _panelmsk = 0x380; // mask for panel field
      constexpr static uint16_t _panelsft = 7; // shift for panel field
      constexpr static uint16_t _planemsk = 0xFC00; // mask for plane field
      constexpr static uint16_t _planesft = 10; // shift for plane field
      constexpr static uint16_t _nstraws = 96; // number of straws per panel
      constexpr static uint16_t _nlayers = 2; // number of layers per panel ; do we need it, see below
      constexpr static uint16_t _npanels = 6; // number of panels per plane
      constexpr static uint16_t _nplanes = 36; // number of planes (36 or 40 for cd3)
      constexpr static uint16_t _invalid = 0xFFFF; // invalid identifier

    public:
      // test values
      static bool validStraw(uint16_t istraw) { return istraw < _nstraws; }
      static bool validPanel(uint16_t ipanel) { return ipanel < _npanels; }
      static bool validPlane(uint16_t iplane) { return iplane < _nplanes; }

      explicit StrawId(std::string const& asstring);

      StrawId(): _sid(_invalid) {}

      // construct from fields
      StrawId( uint16_t plane,
	  uint16_t panel,
	  uint16_t straw);

      // No automatic conversion of uint16_t to StrawId.
      explicit StrawId(uint16_t sid):
        _sid(sid){
        valid();
      }

      // Use compiler-generated copy c'tor, copy assignment, and d'tor.

      // test validity
      bool valid() const { return validStraw(getStraw()) &&
	validPanel(getPanel()) && validPlane(getPlane()); }

      StrawId strawId2() const { return static_cast<StrawId>(_sid); }

      uint16_t asUint16() const { return _sid;}

      uint16_t getPlane() const{
	return (_sid & _planemsk) >> _planesft;
      }

      StrawId getPlaneId() const{
	return static_cast<StrawId>(_sid & _planemsk);
      }

      StrawId getPanelId() const{
	return static_cast<StrawId>(_sid & (_planemsk+_panelmsk) );
      }

      // Id of the first straw in a layer assuming 0 is in 0th etc...
      StrawId getLayerId() const{
	return static_cast<StrawId>( (_sid & (_planemsk+_panelmsk)) +
                                     (_sid & _strawmsk)%_nlayers);
      }

      uint16_t getPanel() const{
	return (_sid & _panelmsk) >> _panelsft;
      }

      uint16_t getStraw() const{
	return (_sid & _strawmsk);
      }

      uint16_t getLayer() const{
	return _sid % _nlayers;
      }

      uint16_t getStation() const{
	return floor(getPlane()/_nlayers);
      }

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
	return _sid % _nlayers;
      }

      uint16_t station() const{
	return floor(plane()/_nlayers);
      }

      bool operator==( StrawId const& rhs) const{
	return ( _sid == rhs._sid );
      }

      bool operator!=( StrawId const& rhs) const{
	return !( *this == rhs);
      }

      bool operator<( StrawId const& rhs) const{
        return ( _sid < rhs._sid);
      }
      bool operator>( StrawId const& rhs) const{
        return ( _sid > rhs._sid);
      }

      bool samePlane( StrawId const& sid) const{
	return ((_sid & _planemsk) == (sid._sid & _planemsk));
      }

      bool samePanel( StrawId const& sid) const{
	return ((_sid & _panelmsk) == (sid._sid & _panelmsk));
      }

      friend std::ostream& operator<<(std::ostream& ost,
                                      const StrawId& s ){
	ost << s.plane() << "_"
            << s.panel() << "_"
            << s.straw();
        return ost;
      }

      // friend std::ostream& operator<<(std::ostream& ost,
      //                                 const StrawId& s ){
      //   ost << "StrawId with Plane " << std::setw(2) << s.plane() << " "
      //       << "Panel " << std::setw(1) << s.panel() << " "
      //       << "Straw " << std::setw(2) << s.straw() << std::endl;
      //   return ost;
      // }

// qualify how close 2 panels are by their Z separation.  This needs to be a logical
// separation, in case there are alignment constants applied
    enum isep{same=0,plane1,station1,station2,station3,apart};

    isep separation(StrawId const& other) const { // was PanelId
      isep retval=apart;
      // same station
      if(other.getPlane()/2 == getPlane()/2){
        int pln1 = getPanel()%2;
        int pln2 = other.getPanel()%2;
        int dp = pln2 - pln1;
        if(other.getPlane() == getPlane()){
          if(dp == 0)
            retval = same;
          else
            retval = plane1;
        } else {
          int dd = other.getPlane() - getPlane();
          if(dp == 0)
            retval = station2;
          else if(dd*dp>0)
            retval = station3;
          else
            retval = station1;
        }
      }
      return retval;
    }

    private:
      // fill fields
      void setStraw(uint16_t istraw);
      void setPanel(uint16_t ipanel);
      void setPlane(uint16_t iplane);

  };
  inline std::ostream& operator<<(std::ostream& ost,
                                  const StrawId& s );

}
#endif /* TrackerGeom_StrawId_hh */
