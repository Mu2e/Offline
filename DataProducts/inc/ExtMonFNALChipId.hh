#ifndef TrackerGeom_ExtMonFNALChipId_hh
#define TrackerGeom_ExtMonFNALChipId_hh

// Identifier of a silicon chip in Mu2e ExtMonFNAL detector.
//
// $Id: ExtMonFNALChipId.hh,v 1.1 2012/08/23 23:41:34 gandr Exp $
// $Author: gandr $
// $Date: 2012/08/23 23:41:34 $
//
// Original author Andrei Gaponenko

#include <ostream>

#include "DataProducts/inc/ExtMonFNALSensorId.hh"

namespace mu2e {

  class ExtMonFNALChipId {
  public:

    ExtMonFNALChipId(const ExtMonFNALSensorId& sensor, int chipCol, int chipRow);

    // Default constructor should not be used by Mu2e code, but it is required by ROOT persistency
    ExtMonFNALChipId() : sensor_(-1), chipCol_(), chipRow_() {}

    const ExtMonFNALSensorId& sensor() const { return sensor_; }
    int   chipCol() const { return chipCol_; }
    int   chipRow() const { return chipRow_; }

    bool operator==( ExtMonFNALChipId const& rhs) const{
      return (sensor_ == rhs.sensor_)&&(chipRow_ == rhs.chipRow_)&&(chipCol_ == rhs.chipCol_);
    }

    bool operator!=( ExtMonFNALChipId const& rhs) const{
      return !(*this == rhs);
    }

    bool operator<( ExtMonFNALChipId const& rhs) const{
      return
        (sensor_ < rhs.sensor_) ||
        ((sensor_ == rhs.sensor_) && ((chipRow_ < rhs.chipRow_) ||
                                      ((chipRow_ == rhs.chipRow_) && (chipCol_ < rhs.chipCol_))
                                      )
         );
    }

  private:
    ExtMonFNALSensorId sensor_;
    int chipCol_;
    int chipRow_;
  };

  std::ostream& operator<<( std::ostream& os, const ExtMonFNALChipId& id);

}
#endif /* TrackerGeom_ExtMonFNALChipId_hh */
