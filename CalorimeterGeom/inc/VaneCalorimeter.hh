#ifndef CalorimeterGeom_VaneCalorimeter_hh
#define CalorimeterGeom_VaneCalorimeter_hh
//
// Hold all geometry and identifier information about
// a Calorimeter.  In order to insulate this class from
// knowledge of databases etc, this class must not know
// how to make itself.
//
// $Id: VaneCalorimeter.hh,v 1.8 2014/08/01 20:57:44 echenard Exp $
// $Author: echenard $
// $Date: 2014/08/01 20:57:44 $
//
// Original author R. Bernstein and Rob Kutschke
//


#include "CalorimeterGeom/inc/BaseCalorimeter.hh"
#include "CalorimeterGeom/inc/Vane.hh"

#include "CLHEP/Vector/ThreeVector.h"

#include <vector>



namespace mu2e {

class VaneCalorimeter: public BaseCalorimeter{

      
      friend class VaneCalorimeterMaker;


      public:

          VaneCalorimeter(){}
          ~VaneCalorimeter(){}


	  //vane components
          int          nVane()            const  {return _nSections; }
	  Vane const&  vane(int i)        const  {return static_cast<Vane const&>(section(i));}

          //few accessors to vane internal composition
	  int          nCrystalX()        const  {return _nCrystalX;}
          int          nCrystalY()        const  {return _nCrystalY;}
          int          crystalX(int crid) const  {return (crid%(_nCrystalX*_nCrystalY))%_nCrystalX;}
          int          crystalY(int crid) const  {return (crid%(_nCrystalX*_nCrystalY))/_nCrystalX;}
          int          vaneId(int crid)   const  {return crid/(_nCrystalX*_nCrystalY);}
          double       innerRadius ()     const  {return _rMin;}
          double       outherRadius()     const  {return _rMax;}
          bool         isTilted()         const  {return _isVaneTilted;}
                  
		  
	
	  //geometry components / print
          virtual bool isInsideCalorimeter(const CLHEP::Hep3Vector &pos)                        const;        
          virtual bool isInsideSection(int iSection, const CLHEP::Hep3Vector &pos)              const;
	  virtual bool isContainedSection(const CLHEP::Hep3Vector &, const CLHEP::Hep3Vector &) const;
	  

	  //crystal id and neighbors component
	  virtual int               crystalIdxFromPosition(const CLHEP::Hep3Vector &pos)  const;
          virtual std::vector<int>  neighborsByLevel(int crystalId, int level)            const; 

	  virtual void print(std::ostream &os = std::cout)                                const;

          


      private:

          int    _nCrystalX;
          int    _nCrystalY;
          double _rMin;
          double _rMax;
          double _shieldHalfThickness;
          double _absorberHalfThickness;
          bool   _isVaneTilted;

   };

}

#endif
