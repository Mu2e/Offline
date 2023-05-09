#ifndef CaloConditions_CalEnergyCalibMaker_hh
#define CaloConditions_CalEnergyCalibMaker_hh

//
// construct a CalEnergyCalib conditions entity
// from fcl or database
// FIXME - currently a place holder
// author: S. Middleton 2022
//
#include "Offline/CaloConditions/inc/CalCalibConstant.hh"
#include "Offline/CaloConfig/inc/CalEnergyCalibConfig.hh" 
#include "Offline/DbTables/inc/CalEnergyCalib.hh"


namespace mu2e {

  class CalEnergyCalibMaker {
    typedef std::shared_ptr<CalCalibConstant> ptr_t;

    public:
      CalEnergyCalibMaker(CalEnergyCalibConfig const& config):_config(config) {};
      ptr_t fromFcl();
      ptr_t fromDb(CalCalibConstant::cptr_t ecalib0);
    
    private:

      const CalEnergyCalibConfig _config;

  };
}


#endif

