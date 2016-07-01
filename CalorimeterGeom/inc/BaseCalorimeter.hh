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

// Mu2e includes
#include "Calorimeter.hh"
#include "BaseCalorimeterInfoGeom.hh"
#include "CaloSection.hh"
#include "Crystal.hh"

//CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"




namespace mu2e {
    

    class BaseCalorimeter: public Calorimeter {

	friend class VaneCalorimeterMaker;
	friend class DiskCalorimeterMaker;


	public:

            BaseCalorimeter();
            virtual ~BaseCalorimeter() {}

	   
	    //calo type
	    virtual const CaloType& caloType() const { return _caloType;}

	    //geomInfo accessor
	    virtual const BaseCalorimeterInfoGeom& caloGeomInfo()  const  {return _caloGeomInfo;}
	    
 	    
	    // Coordinate position and transformation  - The REFERENCE is _Mu2e_ frame
	    // "from" methods: from XXX to Mu2e 
            // "to" methods  : from Mu2e  to XXX     

	    virtual const CLHEP::Hep3Vector&  origin()  const  {return _origin;}
	    virtual CLHEP::Hep3Vector toCrystalFrame(  int crystalId, const CLHEP::Hep3Vector &pos)   const;
	    virtual CLHEP::Hep3Vector toSectionFrame(  int sectionId, const CLHEP::Hep3Vector &pos)   const;	    
	    virtual CLHEP::Hep3Vector toSectionFrameFF(int sectionId, const CLHEP::Hep3Vector &pos)   const;	    
	    virtual CLHEP::Hep3Vector toTrackerFrame(  const CLHEP::Hep3Vector &pos)                  const;

	    virtual CLHEP::Hep3Vector fromCrystalFrame(  int crystalId, const CLHEP::Hep3Vector &pos) const;
	    virtual CLHEP::Hep3Vector fromSectionFrame(  int sectionId, const CLHEP::Hep3Vector &pos) const;
	    virtual CLHEP::Hep3Vector fromSectionFrameFF(int sectionId, const CLHEP::Hep3Vector &pos) const;
	    virtual CLHEP::Hep3Vector fromTrackerFrame(  const CLHEP::Hep3Vector &pos)                const;
            


  	    //calo sections
	    virtual unsigned int              nSection()                                 const  {return _nSections;}
	    virtual const CaloSection&        section(int i)                             const  {return *_sections.at(i);}

            	    
  	    //crystal / readout section
            virtual unsigned int              nCrystal()                                 const  {return _fullCrystalList.size();}
            virtual const Crystal&            crystal(int i)                             const  {return *_fullCrystalList.at(i);}
	            

	    virtual unsigned int              nRO()                                      const  {return _fullCrystalList.size()*_caloGeomInfo.nROPerCrystal();}
	    virtual int                       crystalByRO(int roid)                      const  {return (roid/_caloGeomInfo.nROPerCrystal());}
	    virtual int                       ROBaseByCrystal(int crystalId)             const  {return (crystalId*_caloGeomInfo.nROPerCrystal());}

            
  	    // neighbors, indexing 
	    virtual const std::vector<int>&   neighbors(int crystalId)                   const  {return _fullCrystalList.at(crystalId)->neighbors();}	  
            virtual const std::vector<int>&   nextNeighbors(int crystalId)               const  {return _fullCrystalList.at(crystalId)->nextNeighbors();}	  



	protected:

	    typedef std::shared_ptr<CaloSection> CaloSectionPtr;


	    CaloType                      _caloType;	    
	    unsigned int                  _nSections;
	    std::vector<CaloSectionPtr >  _sections;
	    CLHEP::Hep3Vector             _origin;
	    CLHEP::Hep3Vector             _trackerCenter;

	    BaseCalorimeterInfoGeom       _caloGeomInfo;

	    std::vector<Crystal const*>   _fullCrystalList; //bare pointer is ok, this is just another view for convenience, the objects are managed somewhere else
	    

     };

}    

#endif 
