#ifndef CalorimeterGeom_BaseCalorimeter_hh
#define CalorimeterGeom_BaseCalorimeter_hh
//
// $Id: BaseCalorimeter.hh,v 1.7 2014/02/25 01:09:42 echenard Exp $
// $Author: echenard $
// $Date: 2014/02/25 01:09:42 $
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
#include <boost/shared_ptr.hpp>

// Mu2e includes
#include "Calorimeter.hh"
#include "CaloSection.hh"

//CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"


namespace mu2e {

    

    class BaseCalorimeter: public Calorimeter{

       friend class VaneCalorimeterMaker;
       friend class DiskCalorimeterMaker;
       friend class HybridCalorimeterMaker;
      
       public:

           BaseCalorimeter()  {}
           virtual ~BaseCalorimeter() {}


	   virtual unsigned int             nSection()        const  {return _nSections;}
	           CaloSection       const& section(size_t i) const  {return *_sections.at(i);}
	   virtual CLHEP::Hep3Vector const& origin()          const  {return _origin;}
	   virtual CLHEP::Hep3Vector const& crystalShift()    const  {return _crystalShift;}


           virtual unsigned int nRO()                         const;
           virtual unsigned int nCrystal()                    const;
           virtual unsigned int nROPerCrystal()               const  {return _nROPerCrystal;}

	   virtual int crystalByRO(int roid)                  const  {return (roid/_nROPerCrystal);}
	   virtual int ROBaseByCrystal(int crystalId)         const  {return (crystalId*_nROPerCrystal);}
	           int caloSectionId(int crystalId)           const;
                   int localCrystalId(int crystalId)          const;

	   virtual CLHEP::Hep3Vector crystalOrigin(int crystalId) const;
	   virtual CLHEP::Hep3Vector localCrystalOrigin(int crystalId) const;
           virtual CLHEP::Hep3Vector crystalAxis(int crystalId) const ;
	   virtual CLHEP::Hep3Vector toCrystalFrame(int crystalId, CLHEP::Hep3Vector const& pos) const ;
	   virtual CLHEP::Hep3Vector toSectionFrame(int sectionId, CLHEP::Hep3Vector const& pos) const ;
	   virtual CLHEP::Hep3Vector fromSectionFrame(int sectionId, CLHEP::Hep3Vector const& pos) const ;

	   //a few accessors for convenience
           virtual double wrapperThickness()   const  {return _wrapperThickness;}
           virtual double shellThickness()     const  {return _shellThickness;}
           virtual double caseThickness()      const  {return _caseThickness;}
           virtual double roHalfSize()         const  {return _roHalfTrans;}
           virtual double roHalfThickness()    const  {return _roHalfThickness;}
           virtual double enveloppeInRadius()  const  {return _enveloppeInRadius;}
           virtual double enveloppeOutRadius() const  {return _enveloppeOutRadius;}
           virtual double enveloppeZ0()        const  {return _enveloppeZ0;}
           virtual double enveloppeZ1()        const  {return _enveloppeZ1;}
 
           virtual int    nPipes()             const  {return _nPipes;}
           virtual double pipeRadius()         const  {return _pipeRadius;}
           virtual double pipeThickness()      const  {return _pipeThickness;}
           virtual std::vector<double> const& pipeTorRadius() const  {return _pipeTorRadius;}


	   //keep for backward compatibility, they will go away in the future,
	   virtual double getNonuniformity()  const  {return _nonUniformity; }
	   virtual double getTimeGap()        const  {return _timeGap; }
	   virtual double getElectronEdep()   const  {return _electronEdep; }
	   virtual double getElectronEmin()   const  {return _electronEmin; }
	   virtual double apdMeanNoise()      const  {return _apdMeanNoise;}
	   virtual double apdSigmaNoise()     const  {return _apdSigmaNoise;}
	   virtual double lysoLightYield()    const  {return _lysoLightYield;}
	   virtual double apdQuantumEff()     const  {return _apdQuantumEff;}
	   virtual double apdCollectEff()     const  {return _lightCollectEffAPD;}





       protected:

	   typedef boost::shared_ptr<CaloSection> CaloSectionPtr;

	   unsigned int                 _nSections;
	   std::vector<CaloSectionPtr > _sections;
	   CLHEP::Hep3Vector            _origin;
 	   CLHEP::Hep3Vector            _crystalShift;
	   unsigned int                 _nCrystalTot;

	   unsigned int                 _nROPerCrystal;
           double                       _wrapperThickness;
           double                       _roHalfTrans;
           double                       _roHalfThickness;
           double                       _shellThickness;
           double                       _caseThickness;

           double                       _enveloppeInRadius;
           double                       _enveloppeOutRadius;
           double                       _enveloppeZ0;
           double                       _enveloppeZ1;
	   
	   unsigned int                 _nPipes;
	   double                       _pipeRadius;
	   double                       _pipeThickness;
	   std::vector<double>          _pipeTorRadius;



       private:

          //all these will go away in the future        
          double _nonUniformity;
          double _timeGap;
          double _electronEdep; // energy deposition of charged particle crossing APD
          double _electronEmin; // minimum energy deposition of charged particle crossing APD

          double _apdMeanNoise; //MeV
          double _apdSigmaNoise;//MeV

          double _lysoLightYield;
          double _apdQuantumEff;//quantum efficiency for Hamamatsu S8664-1010 for a radiation wavelenght of 402nm (typical of lyso)
          double _lightCollectEffAPD;//light collection efficiency for 30 × 30 mm2 area of the crystal efficiency for Hamamatsu S8664-1010
     };

}    

#endif /* CalorimeterGeom_BaseCalorimeter_hh */
