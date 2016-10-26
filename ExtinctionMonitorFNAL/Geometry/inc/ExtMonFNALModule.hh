#ifndef EXTMONFNALMODULE_HH
#define EXTMONFNALMODULE_HH

#include <vector>

#include "canvas/Persistency/Common/Wrapper.h"

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALPixelChip.hh"
#include "DataProducts/inc/ExtMonFNALPixelId.hh"
#include "DataProducts/inc/ExtMonFNALModuleId.hh"
#include "CLHEP/Vector/TwoVector.h"

namespace mu2e {

  namespace ExtMonFNAL { class ExtMon; }
  namespace ExtMonFNAL { class ExtMonMaker; }

  class ExtMonFNALModule {
  public:

    std::vector<double> sensorHalfSize() const { return sensorHalfSize_; }
    std::vector<double> chipHalfSize() const { return chipHalfSize_; }

    // Returns default-constructed ExtMonFNALPixelId for out of range (x,y) inputs.
    ExtMonFNALPixelId findPixel(ExtMonFNALModuleId mid, double xModule, double yModule) const;

    // No protection for invalid pixel id
    CLHEP::Hep2Vector moduleCoordinates(const ExtMonFNALPixelId& pix) const;

    const ExtMonFNALPixelChip& chip() const { return chip_; }
    unsigned int nxChips() const;
    unsigned int nyChips() const;

    ExtMonFNALModule(const ExtMonFNALPixelChip& chip, const std::vector<double>& hs)
      : chip_(chip), sensorHalfSize_(hs)
    {}
   // Required by genreflex persistency
    ExtMonFNALModule() {}

    //----------------------------------------------------------------
  private:
    ExtMonFNALPixelChip chip_;
    std::vector<double> chipHalfSize_;
    std::vector<double> sensorHalfSize_;

    template<class T> friend class art::Wrapper;

    friend class ExtMonFNAL::ExtMon;
    friend class ExtMonFNAL::ExtMonMaker;
  };

  //================================================================

} // namespace mu2e

#endif/*EXTMONFNALMODULE_HH*/
