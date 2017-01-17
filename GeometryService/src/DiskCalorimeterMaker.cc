//
// Make a Calorimeter.
//
// original authors  Bertrand Echenarrd

// Disk geometry
//
//  crystals + readout are places in the wrapping material
//  wrapper are placed in the disk
//  disk is placed into the casing, but the casing has only an inside / outside radii thickness, no front back thickness (no cover on the crystals)
//  calibration system (pipes) in front of the crystals
//  pipes + casing form a disk
//
//  cross section in z:      pipes - wrapper - crystal - readout - readout card
//  cross section in radius: casing - wrapper - crystal - wrapper - crystal - ... - wrapper - casing



//  Coordinate systems

//    Intro: We model the crtsyals as polyhedra, so they can be either squares or hexagons. In G4,
//           polyhedra have the origin of the coordinate system at their _base_, and are always
//           oriented along the z direction.
//           On the other hand, the disk have origin of the coordinate system at their _center_
//
//           This choice is the easiest for creation of the geometry and reconstruction.
//           It is also convenient when analyzing track/calo matching to have a coordinate system
//           at the front face of the active volume, the crystals in our case. This one is available
//           in addiiton to the base choice.

//           The idea is to go from the calorimeter -> disk frame -> crystal frame to navigate the geometry


//    Coordinate system
//
//           - crystal:
//                      origin:      base of the crystal
//                      orientation: crystal oriented along the z direction (direction front face - back face)
//                                   the front (back) face is at z=0 (z=crystal length), see above about polyhedra.
//           - disk:
//                      origin:      center of the disk
//                      orientation: along the z axis
//                      other:       extra coordinate system placed at the front of the crystal - used for track-calo matching
//
//           - calorimeter:
//                      origin:      center of front face, point closest to the tracker along z
//
//         Note: the "disk" is defined as whatever volume containing the crystal is rotated in G4. For example, if
//               the crystal + calibration systems are rotated/placed as a single volume, then this is the disk.
//               if the crystals and the calibration systems are rotated/placed separately, then the disk contains only the crystals
//
//         Note: The extra disk coordinate system is suffixed "FF" for FrontFace
//
//    Placement:
//
//         The position of the crystals in the disks (see disk.cc) = the position of the crystal w.r.t disk origin, given by
//         diskOriginToCrystalOrigin:
//              diskOriginToCrystalOrigin(0,0,pipeRadius - crystalHalfLength - roHalfThickness - roElecHalfZ);


// There is a git tag (ef94504f51edbbfeb54a5e63651856bdf5c0a60d) that has the code for a generic placement of the disk origin.
// This is however more complicated, and I think this choice is the best at the moment.



// C++ includes
#include <math.h>
#include <memory>


// Mu2e includes
#include "cetlib/exception.h"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/Disk.hh"
#include "CalorimeterGeom/inc/Crystal.hh"
#include "GeometryService/inc/DiskCalorimeterMaker.hh"

// other includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/RotationX.h"
#include "CLHEP/Vector/RotationY.h"
#include "CLHEP/Vector/RotationZ.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Matrix/DiagMatrix.h"




namespace mu2e {


    DiskCalorimeterMaker::DiskCalorimeterMaker(SimpleConfig const& config, double solenoidOffset)
    {

          _calo = std::unique_ptr<DiskCalorimeter>(new DiskCalorimeter());
          _calo->_caloType = Calorimeter::CaloType::disk;

          _calo->_nSections  = config.getInt("calorimeter.numberOfDisks");
          config.getVectorDouble("calorimeter.diskInnerRadius",  _calo->_diskInnerRadius, _calo->_nSections);
          config.getVectorDouble("calorimeter.diskOuterRadius",  _calo->_diskOuterRadius, _calo->_nSections);
          config.getVectorDouble("calorimeter.diskRotationAngle",_calo->_diskRotAngle,    _calo->_nSections);
          config.getVectorDouble("calorimeter.diskSeparation",   _calo->_diskSeparation,  _calo->_nSections);


          //Fill the Common Calo Data
          _calo->_caloGeomInfo.crystalNedges(      config.getInt(   "calorimeter.crystalNumEdges"));
          _calo->_caloGeomInfo.crystalShift(       config.getBool(  "calorimeter.crystalShift"));
          _calo->_caloGeomInfo.crystalHalfTrans(   config.getDouble("calorimeter.crystalHalfTrans") );
          _calo->_caloGeomInfo.crystalHalfLength(  config.getDouble("calorimeter.crystalHalfLong") );
          _calo->_caloGeomInfo.wrapperThickness(   config.getDouble("calorimeter.crystalWrapperThickness") );
          
	  _calo->_caloGeomInfo.nROPerCrystal(      config.getInt(   "calorimeter.readoutPerCrystal"));
          _calo->_caloGeomInfo.roHalfTrans(        config.getDouble("calorimeter.readoutHalfTrans") );
          _calo->_caloGeomInfo.roHalfThickness(    config.getDouble("calorimeter.readoutHalfThickness") );
          _calo->_caloGeomInfo.roElecHalfX(        config.getDouble("calorimeter.readoutElecHalfX") );
          _calo->_caloGeomInfo.roElecHalfY(        config.getDouble("calorimeter.readoutElecHalfY") );
          _calo->_caloGeomInfo.roElecHalfZ(        config.getDouble("calorimeter.readoutElecHalfZ") );
	  _calo->_caloGeomInfo.crateRadiusIn(      config.getDouble("calorimeter.crateInnerRadius"));
	  _calo->_caloGeomInfo.crateRadiusOut(     config.getDouble("calorimeter.crateOuterRadius"));
 	  _calo->_caloGeomInfo.crateHalfLength(    config.getDouble("calorimeter.crateHalfLength"));
         
	  _calo->_caloGeomInfo.caseThickness(      config.getDouble("calorimeter.caseThickness") );

          _calo->_caloGeomInfo.envelopeInRadius(   config.getDouble("calorimeter.caloMotherInRadius") );
          _calo->_caloGeomInfo.envelopeOutRadius(  config.getDouble("calorimeter.caloMotherOutRadius") );
          _calo->_caloGeomInfo.envelopeZ0(         config.getDouble("calorimeter.caloMotherZ0") );
          _calo->_caloGeomInfo.envelopeZ1(         config.getDouble("calorimeter.caloMotherZ1") );
          _calo->_caloGeomInfo.refractiveIndex(    config.getDouble("calorimeter.refractiveIndex") );
          _calo->_caloGeomInfo.crystalDecayTime(   config.getDouble("calorimeter.crystalDecayTime") );


          _calo->_caloGeomInfo.nPipes(             config.getInt("calorimeter.nPipes",0));
          _calo->_caloGeomInfo.pipeRadius(         config.getDouble("calorimeter.pipeRadius",5) );
          _calo->_caloGeomInfo.pipeThickness(      config.getDouble("calorimeter.pipeThickness",0.5) );

          std::vector<double> temp;
          config.getVectorDouble("calorimeter.pipeTorRadius", temp, _calo->_caloGeomInfo.nPipes());
          _calo->_caloGeomInfo.pipeTorRadius( temp );



          //THE CALORIMETER ORIGIN IS TAKEN AS THE POINT CLOSEST TO THE TRACKER IN MU2E COORDINATES
          double zTrackerCenter  = config.getDouble("mu2e.detectorSystemZ0");
          double xOrigin         = -config.getDouble("mu2e.solenoidOffset");
          double zCaloFront      = config.getDouble("calorimeter.calorimeterZFront");
          double zCaloStart      = config.getDouble("calorimeter.caloMotherZ0");
          double zCaloEnd        = config.getDouble("calorimeter.caloMotherZ1");
          _calo->_origin         = CLHEP::Hep3Vector(xOrigin,0,zCaloFront);
          _calo->_center         = CLHEP::Hep3Vector(xOrigin,0,0.5*(zCaloStart+zCaloEnd));
          
          _calo->_trackerCenter  = CLHEP::Hep3Vector(xOrigin,0,zTrackerCenter);

          
          
          //Get the volume of the solid with hexagonal / rectangular base
          double hl              = _calo->_caloGeomInfo.crystalHalfLength();
          double ht              = _calo->_caloGeomInfo.crystalHalfTrans();
          double cryVolume       = (_calo->_caloGeomInfo.crystalNedges()==4) ? 8*ht*ht*hl : 6.9282032*ht*ht*hl;
          _calo->_caloGeomInfo.crystalVolume(cryVolume);


          _verbosityLevel = config.getInt("calorimeter.verbosityLevel",0);


          //make sure above information is consistent
          CheckIt();

          // Create Disks
          MakeDisks();


    }


    DiskCalorimeterMaker::~DiskCalorimeterMaker() {}




    void DiskCalorimeterMaker::MakeDisks(void)
    {

        int    crystalNedges      = _calo->_caloGeomInfo.crystalNedges();
        bool   crystalShift       = _calo->_caloGeomInfo.crystalShift();
        double crystalHalfLength  = _calo->_caloGeomInfo.crystalHalfLength();
        double crystalHalfTrans   = _calo->_caloGeomInfo.crystalHalfTrans();
        double roHalfThickness    = _calo->_caloGeomInfo.roHalfThickness();
        double roElecHalfZ        = _calo->_caloGeomInfo.roElecHalfZ();
        double caseThickness      = _calo->_caloGeomInfo.caseThickness();
        double wrapperThickness   = _calo->_caloGeomInfo.wrapperThickness();
        double pipeRadius         = _calo->_caloGeomInfo.pipeRadius();
        double crateHalfLength    = _calo->_caloGeomInfo.crateHalfLength();


        double diskHalfZLength    = crystalHalfLength + roHalfThickness + roElecHalfZ + 0.5*wrapperThickness + pipeRadius;
        double crystalCellRadius  = crystalHalfTrans  + wrapperThickness;


        // This is where the offsets between the different coordinate systems are set, see Note for full explanation
        // Seriously, read the note at the top before changing this! Really!
        CLHEP::Hep3Vector diskOriginToCrystalOrigin(0,0,pipeRadius - crystalHalfLength - roHalfThickness - roElecHalfZ + 0.5*wrapperThickness);



        for (unsigned int idisk=0; idisk<_calo->_nSections; ++idisk)
        {
            double dR1    = _calo->_diskInnerRadius[idisk] - caseThickness;
            double dR2    = _calo->_diskOuterRadius[idisk] + caseThickness;
            double dZ     = 2.0*diskHalfZLength;

            CLHEP::Hep3Vector size(dR1,dR2,dZ) ;
            CLHEP::Hep3Vector originLocal(0, 0, diskHalfZLength + _calo->_diskSeparation[idisk]);
            
	    CLHEP::Hep3Vector frontFaceCenter = _calo->origin() + originLocal + diskOriginToCrystalOrigin;
            CLHEP::Hep3Vector backFaceCenter  = frontFaceCenter + CLHEP::Hep3Vector(0,0,2.0*(diskHalfZLength-pipeRadius));

	    std::shared_ptr<Disk> thisDisk( new Disk(idisk,_calo->_diskInnerRadius[idisk], _calo->_diskOuterRadius[idisk], size,
                                                     2.0*crystalCellRadius,crystalNedges, crystalShift, crystalHalfLength, diskOriginToCrystalOrigin) );
            _calo->_sections.push_back(thisDisk);

            thisDisk->setOriginLocal(originLocal );
            thisDisk->setOrigin(_calo->origin() + originLocal );
            thisDisk->setRotation(CLHEP::HepRotation::IDENTITY*CLHEP::HepRotationZ(_calo->_diskRotAngle[idisk]) );
            thisDisk->setFrontFaceCenter(frontFaceCenter);
            thisDisk->setBackFaceCenter(backFaceCenter);
            thisDisk->setCrateDeltaZ(crateHalfLength - diskHalfZLength + pipeRadius);            	    
	    thisDisk->setEnvelopeRad(dR1,dR2);



            //fill the full Crystal List / CaloSectionId (direct access for performance optimization)
            int crystalOffset = _calo->_fullCrystalList.size();
            for (int icry=0;icry<thisDisk->nCrystals();++icry)
            {
                 Crystal& thisCrystal = thisDisk->crystal(icry);
                 _calo->_fullCrystalList.push_back(&thisCrystal);

                 //precompute the neighbors in the global frame
                 thisCrystal.setNeighbors(_calo->neighborsByLevel(icry+crystalOffset,1));
                 thisCrystal.setNextNeighbors(_calo->neighborsByLevel(icry+crystalOffset,2));

                 //pre-compute the crystal position in the mu2e frame (aka global frame), taken from BaseCalorimeter.cc
                 CLHEP::Hep3Vector globalPosition = thisDisk->origin() + thisDisk->inverseRotation()*(thisCrystal.localPosition());
                 thisCrystal.setPosition(globalPosition);
            }


            if (_verbosityLevel) std::cout<<"Constructed Disk "<<thisDisk->id()<<":  Rin="<<thisDisk->innerRadius()<<"  Rout="<<thisDisk->outerRadius()
                                          <<" (X,Y,Z)="<<thisDisk->origin()<<"  local_(X,Y,Z)="<<thisDisk->originLocal()
                                          <<"  with "<<thisDisk->nCrystals()<<" crystals"<<std::endl;

            if (_verbosityLevel > 1) thisDisk->print()                     ;

            if (_verbosityLevel > 2)
            {
                double espace          = thisDisk->estimateEmptySpace();
                double diskVolume    = 3.1415926*(thisDisk->outerRadius()*thisDisk->outerRadius()-thisDisk->innerRadius()*thisDisk->innerRadius());
                double crystalVolume = 3.4641016*crystalCellRadius*crystalCellRadius*thisDisk->nCrystals();

                std::cout<<"Estimated empty space between the disks and the crystals "<<std::endl;
                std::cout<<"Inner edge and crystals = "<<espace<<std::endl;
                std::cout<<"Outer edge and crystals = "<<diskVolume-crystalVolume-espace<<" "<<std::endl;
            }
        }



    }



    void DiskCalorimeterMaker::CheckIt(void)
    {

        int    crystalNedges      = _calo->_caloGeomInfo.crystalNedges();
        int    nROPerCrystal      = _calo->_caloGeomInfo.nROPerCrystal();
        double crystalHalfLength  = _calo->_caloGeomInfo.crystalHalfLength();
        double roHalfThickness    = _calo->_caloGeomInfo.roHalfThickness();
        double roHalfTrans        = _calo->_caloGeomInfo.roHalfTrans();
        double caseThickness      = _calo->_caloGeomInfo.caseThickness();
        double wrapperThickness   = _calo->_caloGeomInfo.wrapperThickness();
        double pipeRadius         = _calo->_caloGeomInfo.pipeRadius();

         if( ! (crystalNedges == 4 || crystalNedges==6 ) )
            {throw cet::exception("DiskCaloGeom") << "calorimeter.crystalNumEdges can only be 4 (squares) or 6 (hexagons)\n";}

        if( ! (nROPerCrystal ==1 || nROPerCrystal ==2 || nROPerCrystal ==4) )
            throw cet::exception("DiskCaloGeom") << "calorimeter.crystalReadoutChannelCount can only be 1,2 or 4.\n";

        // Check size of readouts
        if( nROPerCrystal==1 )
        {
            if (roHalfTrans > _calo->_caloGeomInfo.crystalHalfTrans() )
                throw cet::exception("DiskCaloGeom") << "calorimeter.crystalReadoutHalfTrans > calorimeter.crystalHexsize.\n";

        } else {
            if (roHalfTrans > 0.5*_calo->caloGeomInfo().crystalHalfTrans())
                throw cet::exception("DiskCaloGeom") << "calorimeter.crystalReadoutHalfTrans > 0.5*calorimeter.crystalHalfTrans.\n";
        }



        //check envelope dimensions
        for (unsigned int i=0;i<_calo->_nSections;++i)
        {
          if (_calo->_diskInnerRadius[i] > _calo->_diskOuterRadius[i])
              throw cet::exception("DiskCaloGeom") << "calorimeter.diskInnerRadius > calorimeter.diskOuterRadius for index="<<i<<".\n";

          if ( (_calo->_diskOuterRadius[i] + caseThickness) > _calo->_caloGeomInfo.envelopeOutRadius())
                      throw cet::exception("DiskCaloGeom") << "calorimeter.diskOuterRadius larger than calorimeter mother for index="<<i<<".\n";

          if ( (_calo->_diskInnerRadius[i] - caseThickness) < _calo->_caloGeomInfo.envelopeInRadius())
                      {throw cet::exception("DiskCaloGeom") << "calorimeter.diskInnerRadius smaller than calorimeter mother for index="<<i<<".\n";}

        }


        //check disk length and envelope
        double diskLength   = 2.0*( crystalHalfLength + roHalfThickness + wrapperThickness + caseThickness + pipeRadius);
        double calozEnd     = _calo->_origin.z() + diskLength + _calo->_diskSeparation[_calo->_nSections-1];
        double calozBegin   = _calo->_origin.z();

        if (calozBegin < (_calo->_caloGeomInfo.envelopeZ0()-0.1) || calozBegin > _calo->_caloGeomInfo.envelopeZ1())
            throw cet::exception("DiskCaloGeom") << "calorimeter.calorimeterZFront outside calorimeter mother (need 1mm margin for virtual detectors).\n";

        if (calozEnd   > (_calo->_caloGeomInfo.envelopeZ1()-0.1))
            throw cet::exception("DiskCaloGeom") << "calorimeter z-coordinate extends outside calorimeter mother (need 1mm margin for virtual detectors).\n";



        //look at pipes
        for (int i=0;i<_calo->caloGeomInfo().nPipes();++i)
        {
          if ( (_calo->_caloGeomInfo.pipeTorRadius().at(i)- _calo->_caloGeomInfo.pipeRadius()) <   _calo->_caloGeomInfo.envelopeInRadius())
            throw cet::exception("DiskCaloGeom") << "element "<<i<<" of calorimeter.pipeTorRadius is smaller than disk inner radius\n";

          if ( (_calo->_caloGeomInfo.pipeTorRadius().at(i)+ _calo->_caloGeomInfo.pipeRadius()) >  _calo->_caloGeomInfo.envelopeOutRadius())
            throw cet::exception("DiskCaloGeom") << "element "<<i<<" of calorimeter.pipeTorRadius is larger than disk outer radius\n";
        }

    }


}
