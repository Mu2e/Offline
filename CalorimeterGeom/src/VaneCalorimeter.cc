//
// Geometry and identifier info about the VaneCalorimeter.
//
//
// $Id: VaneCalorimeter.cc,v 1.2 2012/11/17 00:06:25 echenard Exp $
// $Author: echenard $
// $Date: 2012/11/17 00:06:25 $
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
	const Vane& vane = getVane( getCaloSectionId(CrystalId) );
	int ic           = getLocalCrystalId(CrystalId);

	CLHEP::Hep3Vector crysLocalPos = vane.getCrystal(ic).position();

	//if you want to coordinates w.r.t the front face of the crystal, uncomment next two lines
	//CLHEP::Hep3Vector shift(-_crystalHL,0,0);
	//crysLocalPos += shift;

	return (vane.getRotation())*(pos-vane.getOrigin())-crysLocalPos;  
    }

    CLHEP::Hep3Vector VaneCalorimeter::toSectionFrame(int sectionId, CLHEP::Hep3Vector const& pos) const 
    {   
	const Vane& vane = getVane(sectionId);
	return (vane.getRotation())*(pos-vane.getOrigin());
    }

    CLHEP::Hep3Vector VaneCalorimeter::fromSectionFrame(int sectionId, CLHEP::Hep3Vector const& pos) const 
    {   
	const Vane& vane = getVane(sectionId);
	return vane.getOrigin() + vane.getInverseRotation()*pos;
    }

    CLHEP::Hep3Vector VaneCalorimeter::getCrystalAxis(int CrystalId) const 
    {
	const Vane& vane = getVane( getCaloSectionId(CrystalId) );
	CLHEP::Hep3Vector vlocal(1,0,0);
	return vane.getInverseRotation()*vlocal;
    }

    CLHEP::Hep3Vector VaneCalorimeter::getCrystalOrigin(int CrystalId) const 
    {          
       const Vane& vane = getVane( getCaloSectionId(CrystalId) );
       int ic           = getLocalCrystalId(CrystalId);

       CLHEP::Hep3Vector crysLocalPos = vane.getCrystal(ic).position();

       //if you want to coordinates w.r.t the front face of the crystal, uncomment next two lines
       //CLHEP::Hep3Vector shift(-_crystalHL,0,0);
       //crysLocalPos += shift;

       return vane.getOrigin() + vane.getInverseRotation()*crysLocalPos; 
    }
        
    CLHEP::Hep3Vector VaneCalorimeter::getLocalCrystalOrigin(int CrystalId) const 
    {          
       const Vane& vane = getVane( getCaloSectionId(CrystalId) );
       int ic           = getLocalCrystalId(CrystalId);

       CLHEP::Hep3Vector crysLocalPos = vane.getCrystal(ic).position();

       //if you want to coordinates w.r.t the front face of the crystal, uncomment next two lines
       //CLHEP::Hep3Vector shift(-_crystalHL,0,0);
       //crysLocalPos += shift;

       return crysLocalPos; 
    }



    std::vector<int> VaneCalorimeter::getNeighbors(int CrystalId, int level) const 
    {

	int iv = getCaloSectionId(CrystalId);
	int ic = getLocalCrystalId(CrystalId);

	int offset = iv*getVane(0).nCrystals();

	std::vector<int> list = getVane(iv).getNeighbors(ic,level);
	transform(list.begin(), list.end(), list.begin(),bind2nd(std::plus<int>(), offset));  
	return list;
    }



}

/*
    CLHEP::Hep3Vector VaneCalorimeter::toSectionFrame(int crystalId, CLHEP::Hep3Vector const& pos) const 
    {   
	const Vane& vane = getVane( getCaloSectionId(CrystalId) );
	return (vane.getRotation())*(pos-vane.getOrigin());
    }

    CLHEP::Hep3Vector VaneCalorimeter::fromSectionFrame(int crystalId, CLHEP::Hep3Vector const& pos) const 
    {   
	const Vane& vane = getVane( getCaloSectionId(CrystalId) );
	return vane.getOrigin() + vane.getInverseRotation()*pos;
    }
*/
