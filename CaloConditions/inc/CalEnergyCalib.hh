#ifndef CaloConditions_CalEnergyCalib_hh
#define CaloConditions_CalEnergyCalib_hh
//
// CalEnergyCalib collects the net response features of crystal
// used in reconstruction 
// FIXME - this is currently just a place holder
// author: S. Middleton 2022
#include <iostream>
#include <vector>
#include <array>
#include "Offline/CalorimeterGeom/inc/Crystal.hh"
#include "Offline/Mu2eInterfaces/inc/ProditionsEntity.hh"
#include "Offline/DataProducts/inc/CaloId.hh"

namespace mu2e {

  class CalEnergyCalib : public ProditionsEntity {
    public:
    
      struct CalEnergyCorr { //the info here forms the basis of energy fix for each unique calo id
        double _scale;
        double _offset;
        CalEnergyCorr();
        CalEnergyCorr(double scale, double offset) : _scale(scale), _offset(offset) {};
      };
      
      typedef std::shared_ptr<CalEnergyCalib> ptr_t;
      typedef std::shared_ptr<const CalEnergyCalib> cptr_t;
      constexpr static const char* cxname = {"CalEnergyCalib"};
      
      // TODO some construction:
      explicit CalEnergyCalib(uint16_t roid) : ProditionsEntity(cxname), _roid(roid){};
      
      virtual ~CalEnergyCalib(){};

      // TODO here there will be accessors and functions
      uint16_t roid(){ return _roid; }
      std::string algName() { return _algName; }
      
      // TODO add more functionality as required
      /*const CalEnergyCorr&  calibrateEnergy(CaloId& id) const { 
        // TODO - needs to look up the id
        CalEnergyCorr corr; 
        return corr;
      };   
      double getPed(CaloId& id){
        // TODO - needs to look up the id
        return 0.0; 
      };*/
      
    void print(std::ostream& os) const {
      os << "CalEnergyCalib parameters: "  << std::endl;
      os << "  readout ID = " << _roid << " " << std::endl;
      os << "  algorithm used = " << _algName << " " << std::endl;
    }

      template<typename T, size_t SIZE>
      void printArray(std::ostream& os, std::string const& name,
          std::array<T,SIZE> const& a) const{};
          
  private:

      uint16_t _roid;
      std::string _algName;
      
  };
}
#endif

