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
    public:
      constexpr static uint16_t _layermsk = 0x1; // mask for layer field
      constexpr static uint16_t _strawmsk = 0x7F; // mask for straw field
      constexpr static uint16_t _panelmsk = 0x380; // mask for panel field
      constexpr static uint16_t _panelsft = 7; // shift for panel field
      constexpr static uint16_t _planemsk = 0xFC00; // mask for plane field
      constexpr static uint16_t _stationmsk = 0xF800; // mask for station field
      constexpr static uint16_t _stationsft = 11; // shift for station field
      constexpr static uint16_t _planesft = 10; // shift for plane field
      constexpr static uint16_t _invalid = 0xFFFF; // invalid identifier
      constexpr static uint16_t _nstraws = 96; // number of straws per panel
      constexpr static uint16_t _nlayers = 2; // number of layers per panel ; do we need it, see below
      constexpr static uint16_t _npanels = 6; // number of panels per plane
      constexpr static uint16_t _nplanes = 40; // number of planes (36 or 40 for cd3)
      constexpr static uint16_t _nupanels = _npanels * _nplanes; // number of unique panels

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

      // various accessors
      uint16_t asUint16() const { return _sid;}

      uint16_t getPlane() const{
	return plane();
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
	return panel();
      }

      uint16_t getStraw() const{
	return straw();
      }

      uint16_t getLayer() const{
	return layer();
      }

      uint16_t getStation() const{
	return station(); 
      }

      uint16_t plane() const{
	return (_sid & _planemsk) >> _planesft;
      }

      uint16_t panel() const{
	return (_sid & _panelmsk) >> _panelsft;
      }

      uint16_t uniquePanel() const{
	return plane()*_npanels + panel();
      }
      
      uint16_t straw() const{
	return (_sid & _strawmsk);
      }

      uint16_t layer() const{
	return (_sid & _layermsk);
      }

      uint16_t station() const{
	return (_sid & _stationmsk) >> _stationsft;
      }

      // operators

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

      friend std::ostream& operator<<(std::ostream& ost,
                                      const StrawId& s ){
	ost << s.plane() << "_"
            << s.panel() << "_"
            << s.straw();
        return ost;
      }

      // helpers

      bool samePlane( StrawId const& sid) const{
	return ((_sid & _planemsk) == (sid._sid & _planemsk));
      }

      bool samePanel( StrawId const& sid) const{
	return ((_sid & _panelmsk) == (sid._sid & _panelmsk));
      }

// qualify how close 2 panels are by their Z separation.  This needs to be a logical
// separation, in case there are alignment constants applied
      enum isep{same=0,plane1,station1,station2,station3,apart};
      isep separation(StrawId const& other) const;

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
