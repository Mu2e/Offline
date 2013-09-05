//
// Hold information about position of a hexagonal cell
//
// Original author B Echenard - P. Ongmongkolkul
//

//Notes: HexMap tesselates a plane with hexagons, we need a ring
// _posUtil handle number conversion between the map and the disk numbering scheme
//
// CrystalShift indicates the shift of the crystal center from the center of 
// the cell in the calo due to the readout and other things. This is usually disk/vane dependent


// C++ includes
#include <iostream>
#include <cmath>

// Mu2e includes
#include "CalorimeterGeom/inc/CaloSection.hh"
#include "CalorimeterGeom/inc/Barrel.hh"
//#include "CalorimeterGeom/inc/HexMap.hh"
//#include "CalorimeterGeom/inc/DiskCrystalPosUtil.hh"

//CLHEP includes
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"


namespace mu2e {

  Barrel::Barrel(int id, double rin, double rout, double cellSize,int nWheels, int nCrystalWheel, CLHEP::Hep3Vector crystalShift) : 
    CaloSection(id, crystalShift),_radiusIn(rin), _radiusOut(rout),
    _cellSize(cellSize),
    _nWheels(nWheels),_nCrystalWheel(nCrystalWheel)
  { 
    fillCrystals(); 
  }

     

  void Barrel::fillCrystals(void)
  {         

    int nCrystals = _nCrystalWheel*_nWheels;
    double angle = CLHEP::twopi/_nCrystalWheel;
    double stepAngle(0.);
    double x,y,z;
    int wheelIndex;
    double cryR;
    for (int i=0; i< nCrystals; ++i){
      stepAngle = double( i % _nCrystalWheel);
      stepAngle = stepAngle*angle;// + 0.5*angle;
      cryR = (_radiusIn + _radiusOut)/2.0;
      x = cryR*std::sin(stepAngle);
      y = cryR*std::cos(stepAngle);
      wheelIndex = int( i/_nCrystalWheel);//int( (i + 1)/_nCrystalWheel);
      z = (double(wheelIndex)*2.0 + 1.0)*_cellSize;
      CLHEP::Hep3Vector pos(x,y,z);
      //pos += _crystalShift;  //see note 
      _crystalList.push_back( Crystal(i,pos) );
    }

    // now precalculate nearest neighbours list for each crystal
    for (int i=0;i<nCrystals;++i){
      _crystalList[i].setNearestNeighbours(findNeighbors(i,1));
    }

  }



  bool Barrel::isInsideBarrel(CLHEP::Hep3Vector const& pos) const
  {
    double x=pos.x(), y=pos.y(), z=pos.z();
    if( z < 0 || z > (_cellSize*_nWheels)) return false;
    double radius = std::sqrt(x*x + y*y);
    if( radius < _radiusIn || radius > _radiusOut) return false;
    
    return true;
  }

  int Barrel::idxFromPosition(CLHEP::Hep3Vector pos) const 
  {
    double x=pos.x(), y=pos.y(), z=pos.z();
    int nWheels = int(z/_cellSize);
    double id = ( double(_nCrystalWheel) )*nWheels;
    
    double phi_step = CLHEP::twopi/double(_nCrystalWheel);
    double phi = std::atan2(y, x);
    int r = int(phi / phi_step);
    
    id += (double(r) );
    
    return id;    
  }




  std::vector<int> Barrel::neighbors(int crystalId, int level) const
  {
    if (level==1) return _crystalList.at(crystalId).nearestNeighbours();
    return findNeighbors(crystalId,level);
  }
      
  std::vector<int> Barrel::findNeighbors(int crystalId, int level) const
  {
    int Z0 = int(crystalId / _nCrystalWheel);
    int R0 = crystalId % _nCrystalWheel;

    std::vector<int> list;
	  
    //Z0-level -> Z0+level with R=R0+level
    for (int i=-level;i<=level;++i) {
      int Z = Z0+i, R = (R0+level) % _nCrystalWheel;
      if (Z < 0 || Z > (_nWheels-1) ) continue;
      list.push_back(Z*_nCrystalWheel + R);
    }

    //R0+level-1 -> R0-level with Z=Z0+level
    for (int i=level-1;i>=-level;--i) {
      int Z = Z0+level, R = (R0+i) % _nCrystalWheel;
      if ( Z > (_nWheels-1)  ) continue;
      list.push_back(Z*_nCrystalWheel + R); 
    }

    //Z0+level-1 -> Z0-level with R=R0-level
    for (int i=level-1;i>=-level;--i) {
      int Z = Z0+i, R = (R0-level) % _nCrystalWheel;
      if (Z < 0 || Z > (_nWheels-1) ) continue;
      list.push_back(Z*_nCrystalWheel + R); 
    }

    //R0-level+1 -> R0-level-1 with Z=Z0-level
    for (int i=-level+1;i<level;++i) {
      int Z = Z0-level, R = (R0+i) % _nCrystalWheel;
      if ( Z < 0) continue;
      list.push_back(Z*_nCrystalWheel + R);
    }
           
    return list;
	   
  }


  // double Barrel::estimateEmptySpace(void) const
  // {
  //   double sum(0),dx(0.02),dy(0.02);
  //   for (double x=0;x<=1.5*_radiusIn;x+=dx){

  //     double y0 = (x<_radiusIn) ? sqrt(_radiusIn*_radiusIn-x*x) : 0;
  //     for (double y=y0;y<=y0+5*_cellSize;y+=dy){
  // 	int mapIdx = _hexMap.indexFromXY(x/_cellSize,y/_cellSize);
  // 	int iCry   = _posUtil.mapToCrystal(mapIdx);
  // 	if (iCry==-1) sum+=dx*dy;		 
  //     }  
  //   }
  //   return 4.0*sum;
  // }


}

