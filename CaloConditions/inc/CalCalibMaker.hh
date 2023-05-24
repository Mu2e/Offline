#ifndef CaloConditions_CalCalibMaker_hh
#define CaloConditions_CalCalibMaker_hh

//
// construct a CalEnergyCalib conditions entity
// from fcl or database
// FIXME - currently a place holder
// author: S. Middleton 2022
//
#include "Offline/CaloConditions/inc/CalCalib.hh"
#include "Offline/CaloConfig/inc/CalCalibConfig.hh"
#include "Offline/DbTables/inc/CalEnergyCalib.hh"

namespace mu2e {

  class CalCalibMaker {
    typedef std::shared_ptr<CalCalib> ptr_t;

    public:
      CalCalibMaker(CalCalibConfig const& config):_config(config) {};
      ptr_t fromFcl();
      ptr_t fromDb(CalEnergyCalib::cptr_t ecalib0);

    private:

      const CalCalibConfig _config;

  };
}

#endif
