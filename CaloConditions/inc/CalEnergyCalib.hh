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
        double scale;
        double offset;   
      };
      typedef std::shared_ptr<CalEnergyCalib> ptr_t;
      typedef std::shared_ptr<const CalEnergyCalib> cptr_t;
      constexpr static const char* cxname = {"CalEnergyCalib"};
      
      // TODO some construction:
      explicit CalEnergyCalib(double value) : ProditionsEntity(cxname), _value(value){};
      
      virtual ~CalEnergyCalib(){};

      // TODO here there will be accessors and functions
      double value(){ return _value; }

      // TODO Function will run calibration routines
      //const CalEnergyCorr&  calibrateEnergy(CaloId& id) const {};   
      //double getPed(CaloId& id){};
      
      void print(std::ostream& os) const;
      void printVector(std::ostream& os, std::string const& name, 
		      std::vector<double> const& a) const;

      template<typename T, size_t SIZE>
      void printArray(std::ostream& os, std::string const& name,
          std::array<T,SIZE> const& a) const;
          
  private:

      double _value;
      
  };
}
#endif

