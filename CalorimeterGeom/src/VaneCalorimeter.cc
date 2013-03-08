//
// Geometry and identifier info about the VaneCalorimeter.
//
//
// $Id: VaneCalorimeter.cc,v 1.4 2013/03/08 01:22:31 echenard Exp $
// $Author: echenard $
// $Date: 2013/03/08 01:22:31 $
//
// Original author R. Bernstein and Rob Kutschke
//
//C++ includes
#include <algorithm>

//mu2e includes
#include "CalorimeterGeom/inc/VaneCalorimeter.hh"



namespace mu2e {


    CLHEP::Hep3Vector VaneCalorimeter::toCrystalFrame(int CrystalId, CLHEP::Hep3Vector const& pos) const 
    {   
	const Vane& vane0 = vane( caloSectionId(CrystalId) );
	int ic           = localCrystalId(CrystalId);

	CLHEP::Hep3Vector crysLocalPos = vane0.crystal(ic).position();

	//if you want to coordinates w.r.t the front face of the crystal, uncomment next two lines
	CLHEP::Hep3Vector shift(-_crystalHL,0,0);
	crysLocalPos += shift;

	return (vane0.rotation())*(pos-vane0.origin())-crysLocalPos;  
    }

    CLHEP::Hep3Vector VaneCalorimeter::toSectionFrame(int sectionId, CLHEP::Hep3Vector const& pos) const 
    {   
	const Vane& vane0 = vane(sectionId);
	double crystalUnitWidth  = _crystalHW + _wrapperThickness + _shellThickness;
	double crystalUnitLength = _crystalHL + _wrapperThickness;
	
	CLHEP::Hep3Vector vlocal(_roHalfThickness - crystalUnitLength,
                          (_nCrystalR )*crystalUnitWidth,
                          (_nCrystalZ )*crystalUnitWidth );
	return (vane0.rotation())*(pos-vane0.origin()) + vlocal;
    }

    CLHEP::Hep3Vector VaneCalorimeter::fromSectionFrame(int sectionId, CLHEP::Hep3Vector const& pos) const 
    {   
	const Vane& vane0 = vane(sectionId);
	double crystalUnitWidth  = _crystalHW + _wrapperThickness + _shellThickness;
	double crystalUnitLength = _crystalHL + _wrapperThickness;
	
	CLHEP::Hep3Vector vlocal(_roHalfThickness - crystalUnitLength,
                            (_nCrystalR )*crystalUnitWidth,
                            (_nCrystalZ )*crystalUnitWidth );
	return vane0.origin() + vane0.inverseRotation()*(pos-vlocal);
    }

    CLHEP::Hep3Vector VaneCalorimeter::crystalAxis(int CrystalId) const 
    {
	const Vane& vane0 = vane( caloSectionId(CrystalId) );
	CLHEP::Hep3Vector vlocal(1,0,0);
	return vane0.inverseRotation()*vlocal;
    }

    CLHEP::Hep3Vector VaneCalorimeter::crystalOrigin(int CrystalId) const 
    {          
       const Vane& vane0 = vane( caloSectionId(CrystalId) );
       int ic           = localCrystalId(CrystalId);

       CLHEP::Hep3Vector crysLocalPos = vane0.crystal(ic).position();

       //if you want to coordinates w.r.t the front face of the crystal, uncomment next two lines
       CLHEP::Hep3Vector shift(-_crystalHL,0,0);
       crysLocalPos += shift;

       return vane0.origin() + vane0.inverseRotation()*crysLocalPos; 
    }
        
    CLHEP::Hep3Vector VaneCalorimeter::localCrystalOrigin(int CrystalId) const 
    {          
       const Vane& vane0 = vane( caloSectionId(CrystalId) );
       int ic           = localCrystalId(CrystalId);

       CLHEP::Hep3Vector crysLocalPos = vane0.crystal(ic).position();

       //if you want to coordinates w.r.t the front face of the crystal, uncomment next two lines
       CLHEP::Hep3Vector shift(-_crystalHL,0,0);
       crysLocalPos += shift;

       return crysLocalPos; 
    }



    bool VaneCalorimeter::isInsideVane(int ivane, CLHEP::Hep3Vector const& pos) const 
    {   
       CLHEP::Hep3Vector posInSection = toSectionFrame(ivane, pos);
       
       double xlim = _crystalHL + _wrapperThickness + _roHalfThickness + 0.5;
       double ylim = _nCrystalR*(_crystalHW+_wrapperThickness+_shellThickness) + 0.5;   
       double zlim = _nCrystalZ*(_crystalHW+_wrapperThickness+_shellThickness) + 0.5;   
       
       if (posInSection.x() < -xlim || posInSection.x() > xlim ) return false;      
       if (posInSection.y() < -ylim || posInSection.y() > ylim ) return false;      
       if (posInSection.z() < -zlim || posInSection.z() > zlim ) return false;      

       return true;
    }

    int VaneCalorimeter::crystalIdxFromPosition(CLHEP::Hep3Vector const& pos) const 
    {   

       int offset(0);
       for (int ivane=0;ivane<_nVane;++ivane) {
          if ( isInsideVane(ivane,pos) ) {
                CLHEP::Hep3Vector posInSection = toSectionFrame(ivane, pos);
                return offset + vane(ivane).idxFromPosition(posInSection.y(),posInSection.z());
  	  }
          offset +=vane(ivane).nCrystals();
       }       	  
       return -1;
    }



    bool VaneCalorimeter::isInsideCalorimeter(CLHEP::Hep3Vector const& pos) const 
    {   
       for (int ivane=0;ivane<_nVane;++ivane) if (isInsideVane(ivane,pos)) return true;
       return false;    
    }





    std::vector<int> VaneCalorimeter::neighbors(int CrystalId, int level) const 
    {

	int iv = caloSectionId(CrystalId);
	int ic = localCrystalId(CrystalId);

	int offset = iv*vane(0).nCrystals();

	std::vector<int> list = vane(iv).neighbors(ic,level);
	transform(list.begin(), list.end(), list.begin(),bind2nd(std::plus<int>(), offset));  
	return list;
    }



}

/*
    CLHEP::Hep3Vector VaneCalorimeter::toSectionFrame(int crystalId, CLHEP::Hep3Vector const& pos) const 
    {   
	const Vane& vane0 = vane( caloSectionId(CrystalId) );
	return (vane0.rotation())*(pos-vane0.origin());
    }

    CLHEP::Hep3Vector VaneCalorimeter::fromSectionFrame(int crystalId, CLHEP::Hep3Vector const& pos) const 
    {   
	const Vane& vane0 = vane( caloSectionId(CrystalId) );
	return vane0.origin() + vane0.inverseRotation()*pos;
    }
*/
