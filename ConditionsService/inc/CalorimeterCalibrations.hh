#ifndef ConditionsService_CalorimeterCalibrations_hh
#define ConditionsService_CalorimeterCalibrations_hh
//
// Parameters for calorimeter calibrations.
//
// $Id: CalorimeterCalibrations.hh,v 1.1 2013/03/05 20:33:25 aluca Exp $
// $Author: aluca $
// $Date: 2013/03/05 20:33:25 $
//
// Original author
// A. Luca', G. Pezzullo, G. Tassielli, I. Sarra

// C++ includes.
#include <iostream>

// Mu2e includes.
#include "Mu2eInterfaces/inc/ConditionsEntity.hh"
#include "CLHEP/Vector/ThreeVector.h"

//#include "CalorimeterGeom/inc/Calorimeter.hh"


//#include "EventGenerator/inc/GeneratorBase.hh"

namespace mu2e
{
  class SimpleConfig;


  struct CalorimeterCalibrations: virtual public ConditionsEntity{


    CalorimeterCalibrations ( SimpleConfig const& config );

    double LRUpar0( int& crystalId) const{return _LRUpar0;}

    double LRUpar0Err(int& crystalId) const{return _LRUpar0Err;}

    double LINpar0(int& crystalId) const{return _linpar0;}

    double LINpar0Err(int& crystalId) const{return _linpar0Err;}

    double LINpar1(int& crystalId) const{return _linpar1;}

    double LINpar1Err(int& crystalId) const{return _linpar1Err;}

    double LINpar2(int& crystalId) const{return _linpar2;}

    double LINpar2Err(int& crystalId) const{return _linpar2Err;}

    double LINpar3(int& crystalId) const{return _linpar3;}

    double LINpar3Err(int& crystalId) const{return _linpar3Err;}


    double ROpe(int& roId) const{return _ROpe;}
    double ROpeErr(int& roId) const{return _ROpeErr;}

    double ROfano(int& roId) const{return _ROfano;}

    double ROnoise(int& roId) const{return _ROnoise;}
       
    double ADC2MeV(int & roId)      const{return _ADC2MeV;}

  private:

    // We want to discourage multi-phase construction.
    CalorimeterCalibrations ();


    // LRU parameters
    double _LRUpar0;// % / cm
    double _LRUpar0Err;


    //non-linearity parameters
    double _linpar0;
    double _linpar1;
    double _linpar2;
    double _linpar3;
    double _linpar0Err;
    double _linpar1Err;
    double _linpar2Err;
    double _linpar3Err;

    //RO photo-statistic number
    double _ROpe;//p.e. / MeV
    double _ROpeErr;//p.e. / MeV

    double _ROfano;

    //value of the sigma used to do the Gaussian smearing due to the electronic noise
    double _ROnoise;//MeV
    double _ROnoiseSigma;//MeV

    double _ADC2MeV;//conversion factor between ADC count and MeV
  };

  // Shift left (printing) operator.
  inline std::ostream& operator<<(std::ostream& ost,
				  const CalorimeterCalibrations&  ){
    ost << "( CalorimeterCalibrations: to be implemented "
	<< " )";

    return ost;
  }
}

#endif /* ConditionsService_CalorimeterCalibrations_hh */
