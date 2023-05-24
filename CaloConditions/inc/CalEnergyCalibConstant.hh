#ifndef CaloConditions_CalEnergyCalib_hh
#define CaloConditions_CalEnergyCalib_hh
//
// CalEnergyCalibConstant collects the net response features of crystal
// used in reconstruction 
// FIXME - this is currently just a place holder
// author: S. Middleton 2022
#include <iostream>
#include <vector>
#include <array>
#include "Offline/CalorimeterGeom/inc/Crystal.hh"
#include "Offline/Mu2eInterfaces/inc/ProditionsEntity.hh"
#include "Offline/DbTables/inc/CalEnergyCalib.hh"

namespace mu2e {

  class CalEnergyCalibConstant : public ProditionsEntity {
    public:
    
      
      typedef std::shared_ptr<CalEnergyCalibConstant> ptr_t;
      typedef std::shared_ptr<const CalEnergyCalibConstant> cptr_t;
      constexpr static const char* cxname = {"CalEnergyCalibConstant"};
      
      explicit CalEnergyCalibConstant(uint16_t roid) : ProditionsEntity(cxname), _roid(roid){};
      
      virtual ~CalEnergyCalibConstant(){};

      uint16_t roid(){ return _roid; }
      float ADC2MeV() { return _ADC2MeV; }
      int algName() { return _algName; }

      float  GetCalibConstant(uint16_t roid) const {
         // TODO - needs to look up by ID
        return 1.0;
      };
      
      
    void print(std::ostream& os) const {
      os << "CalEnergyCalibConstant parameters: "  << std::endl;
      os << "  readout ID = " << _roid << " " << std::endl;
      os << "  ADC2MeV = " << _ADC2MeV << " " << std::endl;
      os << "  algorithm used = " << _algName << " " << std::endl;
    }

      template<typename T, size_t SIZE>
      void printArray(std::ostream& os, std::string const& name,
          std::array<T,SIZE> const& a) const{};
          
  private:

      uint16_t _roid;
      float _ADC2MeV;
      int _algName;
      
  };
}
#endif

