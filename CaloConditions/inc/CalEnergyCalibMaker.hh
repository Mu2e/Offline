#ifndef CaloConditions_CalEnergyCalibMaker_hh
#define CaloConditions_CalEnergyCalibMaker_hh

//
// construct a CalEnergyCalib conditions entity
// from fcl or database
//

#include "Offline/CaloConditions/inc/CalEnergyCalib.hh"
#include "Offline/CaloConfig/inc/CalEnergyCalibConfig.hh" 
#include "Offline/DbTables/inc/CalSourceCalibTable.hh"


namespace mu2e {

  class CalEnergyCalibMaker {
  typedef std::shared_ptr<CalEnergyCalib> ptr_t;

  public:
    CalEnergyCalibMaker(CalEnergyCalib const& config):_config(config) {}
    ptr_t fromFcl();
    ptr_t fromDb(CalEnergyCalib::cptr_t ecalib0,
		 CalEnergyCalib::cptr_t ecalib1 );
  
  private:

    const CalEnergyCalibConfig _config;

  };
}


#endif

