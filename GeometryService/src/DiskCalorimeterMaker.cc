//
// Make a Disk Calorimeter.
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

#include "cetlib/exception.h"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/Disk.hh"
#include "CalorimeterGeom/inc/Crystal.hh"
#include "GeometryService/inc/DiskCalorimeterMaker.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/RotationX.h"
#include "CLHEP/Vector/RotationY.h"
#include "CLHEP/Vector/RotationZ.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Matrix/DiagMatrix.h"

#include <math.h>
#include <memory>



namespace mu2e {


    DiskCalorimeterMaker::DiskCalorimeterMaker(SimpleConfig const& config, double solenoidOffset)
    {
          
          std::vector<double> temp;
          
          calo_ = std::unique_ptr<DiskCalorimeter>(new DiskCalorimeter());

          //Fill the Common Calo Data
          calo_->nDisks_  = config.getInt("calorimeter.numberOfDisks");
          calo_->caloInfo_.crystalShift(       config.getBool(  "calorimeter.crystalShift"));
          calo_->caloInfo_.crystalHalfTrans(   config.getDouble("calorimeter.crystalHalfTrans") );
          calo_->caloInfo_.crystalHalfLength(  config.getDouble("calorimeter.crystalHalfLong") );
          calo_->caloInfo_.wrapperThickness(   config.getDouble("calorimeter.crystalWrapperThickness") );
          
	  calo_->caloInfo_.nROPerCrystal(      config.getInt(   "calorimeter.readoutPerCrystal"));
          calo_->caloInfo_.roHalfTrans(        config.getDouble("calorimeter.readoutHalfTrans") );
          calo_->caloInfo_.roHalfThickness(    config.getDouble("calorimeter.readoutHalfThickness") );
          calo_->caloInfo_.roElecHalfX(        config.getDouble("calorimeter.readoutElecHalfX") );
          calo_->caloInfo_.roElecHalfY(        config.getDouble("calorimeter.readoutElecHalfY") );
          calo_->caloInfo_.roElecHalfZ(        config.getDouble("calorimeter.readoutElecHalfZ") );
	  calo_->caloInfo_.crateRadiusIn(      config.getDouble("calorimeter.crateInnerRadius"));
	  calo_->caloInfo_.crateRadiusOut(     config.getDouble("calorimeter.crateOuterRadius"));
 	  calo_->caloInfo_.crateHalfLength(    config.getDouble("calorimeter.crateHalfLength"));
         
	  calo_->caloInfo_.caseThickness(      config.getDouble("calorimeter.caseThickness") );

          calo_->caloInfo_.envelopeInRadius(   config.getDouble("calorimeter.caloMotherInRadius") );
          calo_->caloInfo_.envelopeOutRadius(  config.getDouble("calorimeter.caloMotherOutRadius") );
          calo_->caloInfo_.envelopeZ0(         config.getDouble("calorimeter.caloMotherZ0") );
          calo_->caloInfo_.envelopeZ1(         config.getDouble("calorimeter.caloMotherZ1") );
          calo_->caloInfo_.refractiveIndex(    config.getDouble("calorimeter.refractiveIndex") );
          calo_->caloInfo_.crystalDecayTime(   config.getDouble("calorimeter.crystalDecayTime") );


          calo_->caloInfo_.nPipes(             config.getInt("calorimeter.nPipes",0));
          calo_->caloInfo_.pipeRadius(         config.getDouble("calorimeter.pipeRadius",5) );
          calo_->caloInfo_.pipeThickness(      config.getDouble("calorimeter.pipeThickness",0.5) );
          config.getVectorDouble("calorimeter.pipeTorRadius", temp, calo_->caloInfo_.nPipes());
          calo_->caloInfo_.pipeTorRadius( temp );

          temp.clear();
          config.getVectorDouble("calorimeter.diskInnerRadius",  temp, calo_->nDisks_);
          calo_->geomInfo_.diskInnerRadius( temp );
          temp.clear();
          config.getVectorDouble("calorimeter.diskOuterRadius",  temp, calo_->nDisks_);
          calo_->geomInfo_.diskOuterRadius( temp );
          temp.clear();
          config.getVectorDouble("calorimeter.diskRotationAngle",  temp, calo_->nDisks_);
          calo_->geomInfo_.diskRotAngle( temp );
          temp.clear();
          config.getVectorDouble("calorimeter.diskSeparation",  temp, calo_->nDisks_);
          calo_->geomInfo_.diskSeparation( temp );
          



          //THE CALORIMETER ORIGIN IS TAKEN AS THE POINT CLOSEST TO THE TRACKER IN MU2E COORDINATES
          double zTrackerCenter  = config.getDouble("mu2e.detectorSystemZ0");
          double xOrigin         = -config.getDouble("mu2e.solenoidOffset");
          double zCaloFront      = config.getDouble("calorimeter.calorimeterZFront");
          double zCaloStart      = config.getDouble("calorimeter.caloMotherZ0");
          double zCaloEnd        = config.getDouble("calorimeter.caloMotherZ1");
          
          calo_->geomInfo_.origin( CLHEP::Hep3Vector(xOrigin,0,zCaloFront) );
          calo_->geomInfo_.center( CLHEP::Hep3Vector(xOrigin,0,0.5*(zCaloStart+zCaloEnd)) );          
          calo_->geomInfo_.trackerCenter(CLHEP::Hep3Vector(xOrigin,0,zTrackerCenter));

          
          
          //Get the volume of the solid with hexagonal / rectangular base
          double hl              = calo_->caloInfo_.crystalHalfLength();
          double ht              = calo_->caloInfo_.crystalHalfTrans();
          double cryVolume       = 8*ht*ht*hl;
          calo_->caloInfo_.crystalVolume(cryVolume);


          verbosityLevel_ = config.getInt("calorimeter.verbosityLevel",0);


          //make sure above information is consistent
          CheckIt();

          // Create Disks
          MakeDisks();
    }


    DiskCalorimeterMaker::~DiskCalorimeterMaker() {}




    void DiskCalorimeterMaker::MakeDisks(void)
    {

        bool   crystalShift       = calo_->caloInfo_.crystalShift();
        double crystalHalfLength  = calo_->caloInfo_.crystalHalfLength();
        double crystalHalfTrans   = calo_->caloInfo_.crystalHalfTrans();
        double roHalfThickness    = calo_->caloInfo_.roHalfThickness();
        double roElecHalfZ        = calo_->caloInfo_.roElecHalfZ();
        double caseThickness      = calo_->caloInfo_.caseThickness();
        double wrapperThickness   = calo_->caloInfo_.wrapperThickness();
        double pipeRadius         = calo_->caloInfo_.pipeRadius();
        double crateHalfLength    = calo_->caloInfo_.crateHalfLength();


        double diskHalfZLength    = crystalHalfLength + roHalfThickness + roElecHalfZ + 0.5*wrapperThickness + pipeRadius;
        double crystalCellRadius  = crystalHalfTrans  + wrapperThickness;


        // This is where the offsets between the different coordinate systems are set, see Note for full explanation
        // Seriously, read the note at the top before changing this! Really!
        CLHEP::Hep3Vector diskOriginToCrystalOrigin(0,0,pipeRadius - crystalHalfLength - roHalfThickness - roElecHalfZ + 0.5*wrapperThickness);



        for (int idisk=0; idisk<calo_->nDisks_; ++idisk)
        {
            double innerRadius = calo_->geomInfo_.diskInnerRadius().at(idisk);
            double outerRadius = calo_->geomInfo_.diskOuterRadius().at(idisk);
            double separation = calo_->geomInfo_.diskSeparation().at(idisk);
           
            double dR1    = innerRadius - caseThickness;
            double dR2    = calo_->geomInfo_.diskOuterRadius().at(idisk) + caseThickness;
            double dZ     = 2.0*diskHalfZLength;

            CLHEP::Hep3Vector size(dR1,dR2,dZ) ;
            CLHEP::Hep3Vector originLocal(0, 0, diskHalfZLength + separation);
            
	    CLHEP::Hep3Vector frontFaceCenter = calo_->geomInfo_.origin() + originLocal + diskOriginToCrystalOrigin;
            CLHEP::Hep3Vector backFaceCenter  = frontFaceCenter + CLHEP::Hep3Vector(0,0,2.0*(diskHalfZLength-pipeRadius));

	    std::shared_ptr<Disk> thisDisk( new Disk(idisk,innerRadius, outerRadius, 2.0*crystalCellRadius, 
                                                     crystalShift, diskOriginToCrystalOrigin) );
            calo_->disks_.push_back(thisDisk);



            thisDisk->geomInfo().size(size);

            thisDisk->geomInfo().originLocal(originLocal);
            thisDisk->geomInfo().origin(calo_->geomInfo_.origin() + originLocal);
            thisDisk->geomInfo().rotation(CLHEP::HepRotation::IDENTITY*CLHEP::HepRotationZ( calo_->geomInfo().diskRotAngle().at(idisk)) );
            thisDisk->geomInfo().frontFaceCenter(frontFaceCenter);
            thisDisk->geomInfo().backFaceCenter(backFaceCenter);
            thisDisk->geomInfo().crateDeltaZ(crateHalfLength - diskHalfZLength + pipeRadius);            	    
	    thisDisk->geomInfo().envelopeRad(dR1,dR2);



            //fill the full Crystal List / diskId (direct access for performance optimization)
            int crystalOffset = calo_->fullCrystalList_.size();
            for (int icry=0;icry<thisDisk->nCrystals();++icry)
            {
                 Crystal& thisCrystal = thisDisk->crystal(icry);
                 calo_->fullCrystalList_.push_back(&thisCrystal);

                 //precompute the neighbors in the global frame
                 thisCrystal.setNeighbors(calo_->neighborsByLevel(icry+crystalOffset,1));
                 thisCrystal.setNextNeighbors(calo_->neighborsByLevel(icry+crystalOffset,2));

                 //pre-compute the crystal position in the mu2e frame (aka global frame), taken from BaseCalorimeter.cc
                 CLHEP::Hep3Vector globalPosition = thisDisk->geomInfo().origin() + thisDisk->geomInfo().inverseRotation()*(thisCrystal.localPosition());
                 thisCrystal.setPosition(globalPosition);
            }


            if (verbosityLevel_) std::cout<<"Constructed Disk "<<thisDisk->id()<<":  Rin="<<thisDisk->innerRadius()<<"  Rout="<<thisDisk->outerRadius()
                                          <<" (X,Y,Z)="<<thisDisk->geomInfo().origin()<<"  local_(X,Y,Z)="<<thisDisk->geomInfo().originLocal()
                                          <<"  with "<<thisDisk->nCrystals()<<" crystals"<<std::endl;

            if (verbosityLevel_ > 1) thisDisk->print()                     ;

            if (verbosityLevel_ > 2)
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
        int    nROPerCrystal      = calo_->caloInfo_.nROPerCrystal();
        double crystalHalfLength  = calo_->caloInfo_.crystalHalfLength();
        double roHalfThickness    = calo_->caloInfo_.roHalfThickness();
        double roHalfTrans        = calo_->caloInfo_.roHalfTrans();
        double caseThickness      = calo_->caloInfo_.caseThickness();
        double wrapperThickness   = calo_->caloInfo_.wrapperThickness();
        double pipeRadius         = calo_->caloInfo_.pipeRadius();

        if( ! (nROPerCrystal ==1 || nROPerCrystal ==2 || nROPerCrystal ==4) )
            throw cet::exception("DiskCaloGeom") << "calorimeter.crystalReadoutChannelCount can only be 1,2 or 4.\n";

        // Check size of readouts
        if( nROPerCrystal==1 )
        {
            if (roHalfTrans > calo_->caloInfo_.crystalHalfTrans() )
                throw cet::exception("DiskCaloGeom") << "calorimeter.crystalReadoutHalfTrans > calorimeter.crystalHexsize.\n";

        } else {
            if (roHalfTrans > 0.5*calo_->caloInfo_.crystalHalfTrans())
                throw cet::exception("DiskCaloGeom") << "calorimeter.crystalReadoutHalfTrans > 0.5*calorimeter.crystalHalfTrans.\n";
        }



        //check envelope dimensions
        for (int i=0;i<calo_->nDisks_;++i)
        {
            if (calo_->geomInfo_.diskInnerRadius().at(i) > calo_->geomInfo_.diskOuterRadius().at(i))
                throw cet::exception("DiskCaloGeom") << "calorimeter.diskInnerRadius > calorimeter.diskOuterRadius for index="<<i<<".\n";

            if ( (calo_->geomInfo_.diskOuterRadius().at(i) + caseThickness) > calo_->caloInfo_.envelopeOutRadius())
                        throw cet::exception("DiskCaloGeom") << "calorimeter.diskOuterRadius larger than calorimeter mother for index="<<i<<".\n";

            if ( (calo_->geomInfo_.diskInnerRadius().at(i) - caseThickness) < calo_->caloInfo_.envelopeInRadius())
                        {throw cet::exception("DiskCaloGeom") << "calorimeter.diskInnerRadius smaller than calorimeter mother for index="<<i<<".\n";}

        }


        //check disk length and envelope
        double diskLength   = 2.0*( crystalHalfLength + roHalfThickness + wrapperThickness + caseThickness + pipeRadius);
        double calozEnd     = calo_->geomInfo_.origin().z() + diskLength + calo_->geomInfo_.diskSeparation().at(calo_->nDisks_-1);
        double calozBegin   = calo_->geomInfo_.origin().z();

        if (calozBegin < (calo_->caloInfo_.envelopeZ0()-0.1) || calozBegin > calo_->caloInfo_.envelopeZ1())
            throw cet::exception("DiskCaloGeom") << "calorimeter.calorimeterZFront outside calorimeter mother (need 1mm margin for virtual detectors).\n";

        if (calozEnd   > (calo_->caloInfo_.envelopeZ1()-0.1))
            throw cet::exception("DiskCaloGeom") << "calorimeter z-coordinate extends outside calorimeter mother (need 1mm margin for virtual detectors).\n";



        //look at pipes
        for (int i=0;i<calo_->caloInfo().nPipes();++i)
        {
          if ( (calo_->caloInfo_.pipeTorRadius().at(i)- calo_->caloInfo_.pipeRadius()) <   calo_->caloInfo_.envelopeInRadius())
            throw cet::exception("DiskCaloGeom") << "element "<<i<<" of calorimeter.pipeTorRadius is smaller than disk inner radius\n";

          if ( (calo_->caloInfo_.pipeTorRadius().at(i)+ calo_->caloInfo_.pipeRadius()) >  calo_->caloInfo_.envelopeOutRadius())
            throw cet::exception("DiskCaloGeom") << "element "<<i<<" of calorimeter.pipeTorRadius is larger than disk outer radius\n";
        }

    }


}
