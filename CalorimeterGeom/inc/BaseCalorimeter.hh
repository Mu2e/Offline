#ifndef CalorimeterGeom_BaseCalorimeter_hh
#define CalorimeterGeom_BaseCalorimeter_hh
//
// $Id: BaseCalorimeter.hh,v 1.9 2014/08/01 21:54:46 echenard Exp $
// $Author: echenard $
// $Date: 2014/08/01 21:54:46 $
//
// Base class of a cloarimeter. Hold informations about the sections composing 
// the calorimeterand generic algorithms to navigate between the coordinate systems
//
// Original author B. Echenard
//
// Note 1: conversion of crystal <-> readout id
//         readout_id = crystal_id*nRoPerCrystal ... crystal_id*nRoPerCrystal + nRoPerCrystal-1		 


//C++ includes
#include <vector>
#include <memory>
//#include <boost/shared_ptr.hpp>

// Mu2e includes
#include "Calorimeter.hh"
#include "BaseCalorimeterInfoGeom.hh"
#include "CaloSection.hh"
#include "Crystal.hh"

//CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"


namespace mu2e {

    

    class BaseCalorimeter: public Calorimeter{

	friend class VaneCalorimeterMaker;
	friend class DiskCalorimeterMaker;



	public:

            BaseCalorimeter() : _sections(),_caloGeomInfo(),_fullCrystalList(),_crystalSectionId()  {}
            virtual ~BaseCalorimeter() {}


	    virtual BaseCalorimeterInfoGeom const& caloGeomInfo()                const  {return _caloGeomInfo;}
	            CaloSection             const& section(int i)                const  {return *_sections.at(i);}
		    Crystal                 const& crystal(int i)                const  {return *_fullCrystalList.at(i);}
	            
 	    virtual unsigned int nSection()                                  const  {return _nSections;}
            virtual unsigned int nCrystal()                                  const  {return _fullCrystalList.size();}
            virtual unsigned int nRO()                                       const  {return _fullCrystalList.size()*_caloGeomInfo.nROPerCrystal();}

	    virtual          int crystalByRO(int roid)                       const  {return (roid/_caloGeomInfo.nROPerCrystal());}
	    virtual          int ROBaseByCrystal(int crystalId)              const  {return (crystalId*_caloGeomInfo.nROPerCrystal());}
            virtual          int caloSectionId(int i)                        const  {return _crystalSectionId.at(i);}


            virtual std::vector<int> const& neighbors(int crystalId)         const  {return _fullCrystalList.at(crystalId)->neighborsLevel1();}	  
            virtual std::vector<int> const& nextNeighbors(int crystalId)     const  {return _fullCrystalList.at(crystalId)->neighborsLevel2();}	  
            virtual std::vector<int> const& nextNextNeighbors(int crystalId) const  {return _fullCrystalList.at(crystalId)->neighborsLevel3();}	  


	    //geometry and transformations to/from the calorimeter coordinates	   
	    virtual CLHEP::Hep3Vector const& origin() const  {return _origin;}
	    virtual CLHEP::Hep3Vector toCrystalFrame(int crystalId, CLHEP::Hep3Vector const& pos) const ;
	    virtual CLHEP::Hep3Vector fromCrystalFrame(int crystalId, CLHEP::Hep3Vector const& pos) const ;
	    virtual CLHEP::Hep3Vector toSectionFrame(int sectionId, CLHEP::Hep3Vector const& pos) const ;
	    virtual CLHEP::Hep3Vector fromSectionFrame(int sectionId, CLHEP::Hep3Vector const& pos) const ;
            virtual CLHEP::Hep3Vector crystalAxis(int crystalId) const ;
	    virtual CLHEP::Hep3Vector crystalOrigin(int crystalId) const;
	    virtual CLHEP::Hep3Vector crystalOriginInSection(int crystalId) const;


	protected:

	    typedef std::shared_ptr<CaloSection> CaloSectionPtr;
	    unsigned int                 _nSections;
	    std::vector<CaloSectionPtr > _sections;
	    CLHEP::Hep3Vector            _origin;

	    BaseCalorimeterInfoGeom      _caloGeomInfo;

	    std::vector<Crystal const*>  _fullCrystalList;
	    std::vector<int>             _crystalSectionId;
	    



     };

}    

#endif /* CalorimeterGeom_BaseCalorimeter_hh */
