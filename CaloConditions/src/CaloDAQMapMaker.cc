#include "Offline/CaloConditions/inc/CaloDAQMapMaker.hh"
#include "Offline/DataProducts/inc/CaloConst.hh"
#include <vector>


namespace mu2e {
  typedef std::shared_ptr<CaloDAQMap> ptr_t;

  ptr_t CaloDAQMapMaker::fromFcl() {

    // creat this at the beginning since it must be used,
    // partially constructed, to complete the construction
    auto ptr = std::make_shared<CaloDAQMap>(_config.DIRAC2CaloMap(), _config.Calo2DIRACMap());

    return ptr;

  } // end fromFcl

  ptr_t CaloDAQMapMaker::fromDb(CalRoIDMapDIRACToOffline::cptr_t tdtc,
                                CalRoIDMapOfflineToDIRAC::cptr_t tctd ) {

    // initially fill from fcl to get all the constants
    auto ptr = fromFcl();

    // now fill it up ..
    // For calorimeter local array#1: DIRAC2Calo
    // Loops over NumCaloDIRAC*NumChanDIRAC channels
    //
    int NumDIRACTotChannel(CaloConst::_nRawChannel);
    std::vector<uint16_t> dirac2calo(NumDIRACTotChannel);
    for (int i=0;i<NumDIRACTotChannel;i++){
      dirac2calo[i] = tdtc->rowAt(i).caloRoID();
    }
    ptr->setDIRAC2CaloMap(dirac2calo);
    //
    // For Calo crystals to DIRAC: 674*2*2 values
    //
    int NumCaloTotChannel(CaloConst::_nChannel);
    std::vector<uint16_t> calo2dirac(NumCaloTotChannel);

    for (int i=0;i<NumCaloTotChannel;i++){
      calo2dirac[i] = tctd->rowAt(i).diracID();
    }
    ptr->setCalo2DIRACMap(calo2dirac);

    return ptr;

  } // end fromDb

}
