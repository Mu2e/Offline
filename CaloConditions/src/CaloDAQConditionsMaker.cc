#include "CaloConditions/inc/CaloDAQConditionsMaker.hh"
//#include "DataProducts/inc/StrawId.hh"

#include <vector>

using namespace std;

namespace mu2e {
  typedef std::shared_ptr<CaloDAQConditions> ptr_t;

  ptr_t CaloDAQConditionsMaker::fromFcl() {

    // creat this at the beginning since it must be used,
    // partially constructed, to complete the construction
    auto ptr = std::make_shared<CaloDAQConditions>(_config.DIRAC2CaloMap(), _config.Calo2DIRACMap());

    return ptr;

  } // end fromFcl

  ptr_t CaloDAQConditionsMaker::fromDb(
      DIRACtoCalo::cptr_t tdtc,
      CalotoDIRAC::cptr_t tctd ) {

    // initially fill from fcl to get all the constants
    auto ptr = fromFcl();

    // now fill it up .. 
    // For calorimeter local array#1: DIRAC2Calo
    // Loops over 136*200 channels
    //
    int NumCaloDIRAC=136;
    int NumChanDIRAC=20;
    int NumDIRACTotChannel = NumCaloDIRAC*NumChanDIRAC;
    vector<uint16_t> dirac2calo(NumDIRACTotChannel);
    for (int i=0;i<NumDIRACTotChannel;i++){
      dirac2calo[i] = tdtc->rowAt(i).caloRoId();
    }
    ptr->setDIRAC2CaloMap(dirac2calo);
    //
    // For Calo crystals to DIRAC: 674*2*2 values
    //
    int NumDisks=2;
    int NumSensors=2;
    int NumCaloTotChannel =674*NumDisks*NumSensors;
    vector<uint16_t> calo2dirac(NumCaloTotChannel); 

    for (int i=0;i<NumCaloTotChannel;i++){
      calo2dirac[i] = tctd->rowAt(i).dirac();
    }
    ptr->setCalo2DIRACMap(calo2dirac);
    
    return ptr;

  } // end fromDb

}
