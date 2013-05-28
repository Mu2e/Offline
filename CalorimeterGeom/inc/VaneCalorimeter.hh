#ifndef CalorimeterGeom_VaneCalorimeter_hh
#define CalorimeterGeom_VaneCalorimeter_hh
//
// Hold all geometry and identifier information about
// a Calorimeter.  In order to insulate this class from
// knowledge of databases etc, this class must not know
// how to make itself.
//
// $Id: VaneCalorimeter.hh,v 1.7 2013/05/28 22:11:24 echenard Exp $
// $Author: echenard $
// $Date: 2013/05/28 22:11:24 $
//
// Original author R. Bernstein and Rob Kutschke
//

//C++ includes
#include <vector>
#include <boost/shared_ptr.hpp>

// Mu2e includes
#include "CalorimeterGeom/inc/BaseCalorimeter.hh"
#include "CalorimeterGeom/inc/Vane.hh"

//CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"


using boost::static_pointer_cast;
using boost::shared_ptr;



namespace mu2e {

class VaneCalorimeter: public BaseCalorimeter{

      
      friend class VaneCalorimeterMaker;


      public:

          VaneCalorimeter(){}
          ~VaneCalorimeter(){}


          int nVane()                           const  {return _nSections; }
	  Vane const&  vane(int i)              const  {return static_cast<Vane const&>(section(i));}

                  int    nCrystalR()            const  {return _nCrystalR; }
                  int    nCrystalZ()            const  {return _nCrystalZ; }
          virtual double crystalHalfTrans()     const  {return _crystalHW; }
          virtual double crystalHalfLength()    const  {return _crystalHL; }
          virtual double crystalVolume()        const  {return 8*_crystalHW*_crystalHW*_crystalHL;}

                  bool             isInsideVane(int ivane, CLHEP::Hep3Vector const& pos) const ;
          virtual bool             isInsideCalorimeter(CLHEP::Hep3Vector const& pos) const ;        
	  virtual int              crystalIdxFromPosition(CLHEP::Hep3Vector const& pos) const ;
          virtual std::vector<int> neighbors(int crystalId, int level=1) const;
          virtual double           crystalLongPos(int crystalId, CLHEP::Hep3Vector const& pos) const; 

          double innerRadius ()                 const  {return _rMin;}
          double outherRadius()                 const  {return _rMax;}


        


//keep only for backward compatibility, will disappear in the future.

int crystalByRO(int roid) const              {return (roid/_nROPerCrystal); }
int ROBaseByCrystal(int crystalId) const     {return (crystalId*_nROPerCrystal);}

int crystalVaneByRO(int roid) const { return (roid/_nROPerCrystal)%(_nCrystalZ*_nCrystalR);}
int crystalRByRO(int roid) const {return ((roid/_nROPerCrystal)%(_nCrystalZ*_nCrystalR))/_nCrystalZ;}
int crystalZByRO(int roid) const {return ((roid/_nROPerCrystal)%(_nCrystalZ*_nCrystalR))%_nCrystalZ;}

CLHEP::Hep3Vector crystalOriginByRO(int roid) const { return crystalOrigin(crystalByRO(roid));  }
int vaneByRO(int roid) const {return roid/(_nCrystalZ*_nCrystalR*_nROPerCrystal);}
CLHEP::Hep3Vector crystalAxisByRO(int roid) const { return crystalAxis(crystalByRO(roid));  }

CLHEP::Hep3Vector toVaneFrame(int vaneId, CLHEP::Hep3Vector const& pos) const   {return toSectionFrame(vaneId,pos);}
CLHEP::Hep3Vector fromVaneFrame(int vaneId, CLHEP::Hep3Vector const& pos) const {return fromSectionFrame(vaneId,pos);}




      private:

          int    _nCrystalZ;
          int    _nCrystalR;
          double _rMin;
          double _rMax;
          double _crystalHL;
          double _crystalHW;        
          double _shieldHalfThickness;
          double _absorberHalfThickness;



   };

}

#endif /* CalorimeterGeom_VaneCalorimeter_hh */
