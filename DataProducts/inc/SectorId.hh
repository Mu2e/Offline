#ifndef TrackerGeom_SectorId_hh
#define TrackerGeom_SectorId_hh

//
// Identifier for a sector.
//

//
// $Id: SectorId.hh,v 1.9 2013/03/08 04:31:45 brownd Exp $
// $Author: brownd $
// $Date: 2013/03/08 04:31:45 $
//
// Original author Rob Kutschke
//

#include <ostream>
#include "DataProducts/inc/DeviceId.hh"

namespace mu2e {

  struct SectorId;
  inline std::ostream& operator<<(std::ostream& ost,
                                  const SectorId& s );

  class SectorId{

  public:
// qualify how close 2 sectors are by their Z separation.  This needs to be a logical
// separation, in case there are alignment constants applied
    enum isep{same=0,device,station1,station2,station3,apart};

    SectorId():
      _did(-1),
      _sector(-1){
    }

    SectorId( DeviceId device,
              int sector
              ):
      _did(device),
      _sector(sector){
    }

    // Use compiler-generated copy c'tor, copy assignment, and d'tor.

    const DeviceId& getDeviceId() const {
      return _did;
    }

    int getDevice() const {
      return _did;
    }

    int getSector() const {
      return _sector;
    }

    bool operator==(SectorId const& rhs) const{
      return ( _did == rhs._did && _sector == rhs._sector );
    }

    bool operator!=(SectorId const& rhs) const{
      return !( *this == rhs);
    }

    bool operator < (SectorId const& rhs) const {
      return _did < rhs._did || (_did == rhs._did && _sector < rhs._sector);
    }

    isep separation(SectorId const& other) const;

    friend std::ostream& operator<<(std::ostream& ost,
                                    const SectorId& s ){
      ost << s._did << " " << s._sector;
      return ost;
    }

  private:

    DeviceId _did;
    int      _sector;

  };

}  //namespace mu2e

#endif /* TrackerGeom_SectorId_hh */
