#ifndef CaloConditions_AlignedCalorimeterMaker_hh
#define CaloConditions_AlignedCalorimeterMaker_hh
//
// Make AlignedCalorimeter from fcl or (eventually) database
//
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/CaloConfig/inc/AlignedCalConfig.hh"
#include "Offline/DbTables/inc/CalAlignElement.hh"
#include "Offline/DbTables/inc/CalAlignParams.hh"

namespace mu2e {

  class AlignedCalorimeterMaker {
      typedef std::shared_ptr<DiskCalorimeter> ptr_t;

    public:
      AlignedCalorimeterMaker(const AlignedCalConfig& config):_config(config) {}
      void alignCalorimeter(ptr_t ptr,
                            const std::vector<CalAlignParams>& disk_align_params);
      ptr_t fromFcl();
      ptr_t fromDb(CalAlignDisk::cptr_t cad_p);

    private:
      std::vector<CalAlignParams> readFile(const std::string& fileName, uint16_t nRowMax);
      std::vector<CalAlignParams> readDb(CalAlignElement::cptr_t cptr, uint16_t nRowMax);

      // this object needs to be thread safe,
      // _config should only be initialized once
      const AlignedCalConfig _config;
  };
}

#endif
