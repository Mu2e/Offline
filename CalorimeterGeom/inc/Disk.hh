#ifndef CalorimeterGeom_CaloSection_hh
#define CalorimeterGeom_CaloSection_hh
//
// Hold information about a disk in the calorimter.
//
// Original author B Echenard 
//


#include "CalorimeterGeom/inc/DiskGeomInfo.hh"
#include "CalorimeterGeom/inc/Crystal.hh"
#include "CalorimeterGeom/inc/CrystalMapper.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/Rotation.h"
#include <vector>
#include <memory>


namespace mu2e {

   class Disk {

       public:
           
           Disk(int id, double rin, double rout,double cellSize, bool shiftCrystal,
		const CLHEP::Hep3Vector& diskOriginToCrystalOrigin);
           ~Disk(){};     


	   int                        id()                     const {return id_;}
           
	   int                        nCrystals()              const {return crystalList_.size();}
           const Crystal&             crystal(int i)           const {return crystalList_.at(i);}
                 Crystal&             crystal(int i)                 {return crystalList_.at(i);}   
           
           const DiskGeomInfo&        geomInfo()               const {return geomInfo_;}
                 DiskGeomInfo&        geomInfo()                     {return geomInfo_;}           
	   
           double                     innerRadius()            const {return radiusIn_;}
           double                     outerRadius()            const {return radiusOut_;}
 	   
	   int                        idxFromPosition(double x, double y) const;           
	   std::vector<int>           findLocalNeighbors(int crystalId, int level,bool raw=false) const;            
           std::vector<int>           nearestIdxFromPosition(double x, double y) const; 
           
           double                     estimateEmptySpace() const;
           void                       print(std::ostream& os = std::cout) const;


       private:
       
           void                            fillCrystals(const CLHEP::Hep3Vector &crystalOriginInDisk);
           bool                            isInsideDisk(double x, double y) const;
	   double                          calcDistToSide(const CLHEP::Hep2Vector& a, const CLHEP::Hep2Vector& b) const;
           
           std::vector<Crystal>            crystalList_;
	   int                             id_;
           DiskGeomInfo                    geomInfo_;
	   double                          radiusIn_;
	   double                          radiusOut_;
	   double                          cellSize_;	 
       	   std::shared_ptr<CrystalMapper>  crystalMap_;
           std::vector<int>                mapToCrystal_;
	   std::vector<int>                crystalToMap_;

   };

}

#endif
