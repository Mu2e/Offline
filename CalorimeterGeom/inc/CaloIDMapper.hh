//
// Contains code to convert SiPM ID into crystal ID and vice-versa
//
// conversion of crystal <-> readout id
// readout_id = crystal_id*nSiPMPerCrystal ... crystal_id*nSiPMPerCrystal + nSiPMPerCrystal-1
//
#ifndef CalorimeterGeom_CaloIDMapper_hh
#define CalorimeterGeom_CaloIDMapper_hh

#include "cetlib_except/exception.h"

#include <vector>
#include <map>
#include <string>


namespace mu2e {

    class CaloIDMapper {

       public:
           CaloIDMapper() : nSiPMPerCrystal_(1) {}

           void nSiPMPerCrystal(int value)               {nSiPMPerCrystal_ = value;}
             int  nSiPMPerCrystal()                  const {return nSiPMPerCrystal_;}
           int  crystalIDFromSiPMID(int roid)      const {return roid/nSiPMPerCrystal_;}
           int  SiPMIDFromCrystalID(int crystalId) const {return crystalId*nSiPMPerCrystal_;}
           int  SiPMIdx(int SiPMID)                const {return SiPMID%nSiPMPerCrystal_;}

       private:
          int nSiPMPerCrystal_;
     };

}

#endif
