#ifndef CalorimeterGeom_HybridCalorimeter_hh
#define CalorimeterGeom_HybridCalorimeter_hh
//
// $Id: HybridCalorimeter.hh,v 1.1 2013/09/05 17:11:28 gianipez Exp $
// $Author: gianipez $
// $Date: 2013/09/05 17:11:28 $
//
// Look at Mu2eG4/inc/constructHybridCalorimeter.cc 
// for definition of geometry

// Original author Gianipez


//C++ includes
#include <vector>
#include <boost/shared_ptr.hpp>

// Mu2e includes
#include "CalorimeterGeom/inc/BaseCalorimeter.hh"
#include "CalorimeterGeom/inc/Disk.hh"
#include "CalorimeterGeom/inc/Barrel.hh"

//CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {


  class HybridCalorimeter: public BaseCalorimeter{

       
    friend class HybridCalorimeterMaker;


  public:


    HybridCalorimeter()  {}
    ~HybridCalorimeter() {}

    unsigned int nSections()                const  {return _nSections;}
    Disk const&  disk()                 const  {return static_cast<Disk const&>(section(0));}
    Barrel const&  barrel()             const  {return static_cast<Barrel const&>(section(1));}
    
    virtual double hexCrystalHalfTrans()   const  {return _hexCrystalHalfTrans;}
    virtual double hexCrystalHalfLength()  const  {return _hexCrystalHalfLength;}
    virtual double hexCrystalVolume()      const  {return _hexCrystalVolume;}//
    
    virtual double crystalHalfTrans()      const  {return _crystalHalfTrans;}
    virtual double crystalHalfLength()     const  {return _crystalHalfLength;}
    virtual double barrelCrystalVolume()   const  {return _barrelCrystalVolume;}
    virtual double crystalVolume()         const  {return (_hexCrystalVolume + _barrelCrystalVolume);}

    bool                     isInsideDisk(CLHEP::Hep3Vector const& pos) const ;       	 
    bool                     isInsideBarrel(CLHEP::Hep3Vector const& pos) const ;       	    
    virtual bool             isInsideCalorimeter(CLHEP::Hep3Vector const& pos) const ;       	 	 
    virtual bool             isInsideSection(int iSec, CLHEP::Hep3Vector const& pos) const;
    virtual int              crystalIdxFromPosition(CLHEP::Hep3Vector const& pos) const ;
    virtual std::vector<int> neighbors(int crystalId, int level=1) const; 
    virtual double           crystalLongPos(int crystalId, CLHEP::Hep3Vector const& pos) const; 



  private:

    double _diskInnerRadius;
    double _diskOuterRadius;
    double _diskRotAngle;  
    double _hexCrystalVolume;
    
    double _nWheels;
    double _nCrystalWheel;
    double _barrelInnerRadius;
    double _barrelOuterRadius;
    double _diskToBarrelSeparation;
    double _barrelRotAngle;
    double _barrelCrystalVolume;
    
    double _crystalHalfTrans;
    double _crystalHalfLength;

    double _hexCrystalHalfTrans;
    double _hexCrystalHalfLength;

  };

}    

#endif /* CalorimeterGeom_HybridCalorimeter_hh */
