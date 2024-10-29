#ifndef CaloConditions_CalCalib_hh
#define CaloConditions_CalCalib_hh
//
// CalCalib collects the net response features of crystal
// used in reconstruction
// FIXME - this is currently just a place holder
// author: S. Middleton 2022
#include <iostream>
#include <vector>
#include <array>
#include "Offline/CalorimeterGeom/inc/Crystal.hh"
#include "Offline/Mu2eInterfaces/inc/ProditionsEntity.hh"
#include "Offline/DbTables/inc/CalEnergyCalib.hh"
#include "Offline/CaloConditions/inc/CalCalibPar.hh"

namespace mu2e {

  class CalCalib : public ProditionsEntity {
    public:

      typedef std::shared_ptr<CalCalib> ptr_t;
      typedef std::shared_ptr<const CalCalib> cptr_t;
      typedef std::vector<CalCalibPar> CalibVec;
      constexpr static const char* cxname = {"CalCalib"};

      explicit CalCalib(const CalibVec& cvec) : ProditionsEntity(cxname), _cvec(cvec){};

      virtual ~CalCalib(){};

      const CalCalibPar& calib(std::uint16_t roid) const {
        return _cvec.at(roid);
      }

      float ADC2MeV(std::uint16_t roid) const {
        return _cvec.at(roid).ADC2MeV();
      }

      float timeoffset(std::uint16_t roid) const {
        return _cvec.at(roid).timeOffset();
      }

      void print( std::ostream& ) const;

  private:

      CalibVec _cvec;

  };
}
#endif
