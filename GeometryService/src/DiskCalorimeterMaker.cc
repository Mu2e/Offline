//
// Make a Disk Calorimeter.
//
// original authors  Bertrand Echenarrd

// Disk geometry
//
//  see Mu2eG4/src/constructDiskCalorimeter.cc for the geometry
//    the disk are split in three pieces: front palte, disk case and back plate
//  
//  front plate has calibration pipes
//  disk case has crystal enveloped in wrapper (except the back of the crystal!)
//  back plate has a solid plate with holes and readouts inside + FEE box outside


//    Coordinate system
//
//           - crystal:
//                      origin:      base of the crystal
//                      orientation: crystal oriented along the z direction (direction front face - back face)
//                                   the front (back) face is at z=0 (z=crystal length).
//           - disk:
//                      origin:      center of the disk
//                      orientation: along the z axis
//                      other:       extra coordinate system placed at the front of the crystal - used for track-calo matching
//
//           - calorimeter:
//                      origin:      center of mother volume
//
//         Note: The extra disk coordinate system is suffixed "FF" for FrontFace. The FF z origin is at the same position of the crystal z origin.
//


// For reference, git tag (ef94504f51edbbfeb54a5e63651856bdf5c0a60d)  has generic placement of the disk origin.

#include "cetlib_except/exception.h"
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
          verbosityLevel_ = config.getInt("calorimeter.verbosityLevel",0);
          calo_->caloInfo_.envelopeRadiusIn(   config.getDouble("calorimeter.caloMotherRadiusIn") );
          calo_->caloInfo_.envelopeRadiusOut(  config.getDouble("calorimeter.caloMotherRadiusOut") );
          calo_->caloInfo_.envelopeZ0(         config.getDouble("calorimeter.caloMotherZ0") );
          calo_->caloInfo_.envelopeZ1(         config.getDouble("calorimeter.caloMotherZ1") );

          calo_->nDisks_  = config.getInt("calorimeter.numberOfDisks");
          calo_->caloInfo_.crystalShift(       config.getBool(  "calorimeter.crystalShift"));
          calo_->caloInfo_.crystalXYLength(    config.getDouble("calorimeter.crystalXYLength") );
          calo_->caloInfo_.crystalZLength(     config.getDouble("calorimeter.crystalZLength") );
          calo_->caloInfo_.wrapperThickness(   config.getDouble("calorimeter.wrapperThickness") );
          calo_->caloInfo_.refractiveIndex(    config.getDouble("calorimeter.refractiveIndex") );
          calo_->caloInfo_.crystalDecayTime(   config.getDouble("calorimeter.crystalDecayTime") );
          
	  calo_->caloInfo_.nROPerCrystal(      config.getInt(   "calorimeter.readoutPerCrystal"));
          calo_->caloInfo_.roXLength(          config.getDouble("calorimeter.readoutXLength") );
          calo_->caloInfo_.roYLength(          config.getDouble("calorimeter.readoutYLength") );
          calo_->caloInfo_.roZLength(          config.getDouble("calorimeter.readoutZLength") );
          calo_->caloInfo_.FEEXLength(         config.getDouble("calorimeter.FEEXLength") );
          calo_->caloInfo_.FEEYLength(         config.getDouble("calorimeter.FEEYLength") );
          calo_->caloInfo_.FEEZLength(         config.getDouble("calorimeter.FEEZLength") );
          calo_->caloInfo_.FEEBoxThickness(    config.getDouble("calorimeter.FEEBoxThickness") );
          calo_->caloInfo_.BPHoleXLength(      config.getDouble("calorimeter.BPHoleXLength") );
          calo_->caloInfo_.BPHoleYLength(      config.getDouble("calorimeter.BPHoleYLength") );
          calo_->caloInfo_.BPHoleZLength(      config.getDouble("calorimeter.BPHoleZLength") );
          calo_->caloInfo_.stripThickness(     config.getDouble("calorimeter.stripThickness") );
          calo_->caloInfo_.stripYLength(       config.getDouble("calorimeter.stripYLength") );
          calo_->caloInfo_.coolBPPipeRadius(   config.getDouble("calorimeter.coolBPPipeRadius") );

	  calo_->caloInfo_.caseThicknessIn(      config.getDouble("calorimeter.caseThicknessIn") );
	  calo_->caloInfo_.caseThicknessOut(     config.getDouble("calorimeter.caseThicknessOut") );
	  calo_->caloInfo_.outerRingEdgeZLength( config.getDouble("calorimeter.outerRingEdgeZLength") );
	  calo_->caloInfo_.outerRingEdgeRLength( config.getDouble("calorimeter.outerRingEdgeRLength") );

          calo_->caloInfo_.FPstyrofoamZLength( config.getDouble("calorimeter.FPStyrofoamZLength") );
          calo_->caloInfo_.FPCarbonZLength(    config.getDouble("calorimeter.FPCarbonZLength") );
          calo_->caloInfo_.coolFPPipeRadius(   config.getDouble("calorimeter.coolFPPipeRadius") );
          calo_->caloInfo_.nPipes(             config.getInt("calorimeter.nPipes",0));
          calo_->caloInfo_.pipeRadius(         config.getDouble("calorimeter.pipeRadius",5) );
          calo_->caloInfo_.pipeThickness(      config.getDouble("calorimeter.pipeThickness",0.5) );
          config.getVectorDouble("calorimeter.pipeTorRadius", temp, calo_->caloInfo_.nPipes());
          calo_->caloInfo_.pipeTorRadius( temp );
          calo_->caloInfo_.pipeSeparation(     config.getDouble("calorimeter.pipeSeparation") );
          
	  calo_->caloInfo_.nBoard(                 config.getInt("calorimeter.numberOfBoards") );
	  calo_->caloInfo_.nCrate(                 config.getInt("calorimeter.numberOfCrates") );
	  calo_->caloInfo_.nCrateBeforeSpace(      config.getInt("calorimeter.nCrateBeforeSpace") );
	  calo_->caloInfo_.crateXLength(           config.getDouble("calorimeter.crateXLength"));
	  calo_->caloInfo_.crateYLength(           config.getDouble("calorimeter.crateYLength"));
	  calo_->caloInfo_.crateZLength(           config.getDouble("calorimeter.crateZLength"));
          calo_->caloInfo_.crateFShieldThickness(  config.getDouble("calorimeter.crateFShieldThickness") );
          calo_->caloInfo_.crateBShieldThickness(  config.getDouble("calorimeter.crateBShieldThickness") );
          calo_->caloInfo_.crateTThickness(        config.getDouble("calorimeter.crateTThickness") );
          calo_->caloInfo_.crateSThickness(        config.getDouble("calorimeter.crateSThickness") );
          calo_->caloInfo_.crateFShieldYLength(    config.getDouble("calorimeter.crateFShieldYLength") );
	  calo_->caloInfo_.crateFShieldDeltaZ(     config.getDouble("calorimeter.crateFShieldDeltaZ"));
	  calo_->caloInfo_.crateRadiusIn(          config.getDouble("calorimeter.crateInnerRadius"));
          calo_->caloInfo_.cratephi0(              config.getDouble("calorimeter.cratephi0") );
          calo_->caloInfo_.crateDeltaPhi(          config.getDouble("calorimeter.crateDeltaPhi") );
	  calo_->caloInfo_.radiatorThickness(      config.getDouble("calorimeter.radiatorThickness") );
	  calo_->caloInfo_.activeStripThickness(   config.getDouble("calorimeter.activeStripThickness") );
	  calo_->caloInfo_.passiveStripThickness(  config.getDouble("calorimeter.passiveStripThickness") );
	  

          temp.clear();
          config.getVectorDouble("calorimeter.diskCaseRadiusIn",  temp, calo_->nDisks_);
          calo_->geomInfo_.diskCaseRadiusIn( temp );
          temp.clear();
          config.getVectorDouble("calorimeter.diskCaseRadiusOut",  temp, calo_->nDisks_);
          calo_->geomInfo_.diskCaseRadiusOut( temp );
          temp.clear();
          config.getVectorDouble("calorimeter.diskRotationAngle",temp, calo_->nDisks_);
          calo_->geomInfo_.diskRotAngle( temp );
          temp.clear();
          config.getVectorDouble("calorimeter.diskZ0MotherShift", temp, calo_->nDisks_);
          calo_->geomInfo_.diskZMotherShift( temp );
          
          if (!calo_->caloInfo_.nPipes()) calo_->caloInfo_.FPstyrofoamZLength(0);
          if (!calo_->caloInfo_.nPipes()) calo_->caloInfo_.FPCarbonZLength(0);

          if (!calo_->caloInfo_.nROPerCrystal())
          {
              calo_->caloInfo_.roZLength(0);
              calo_->caloInfo_.FEEZLength(0);
              calo_->caloInfo_.FEEBoxThickness(0);
              calo_->caloInfo_.BPHoleZLength(0);
              calo_->caloInfo_.stripThickness(0);
          }


          //CALORIMETER ORIGIN AND FRONT FACE (FF)
          double zTrackerCenter  = config.getDouble("mu2e.detectorSystemZ0");
          double xOrigin         = -config.getDouble("mu2e.solenoidOffset");
          double zCaloStart      = config.getDouble("calorimeter.caloMotherZ0");
          double zCaloEnd        = config.getDouble("calorimeter.caloMotherZ1");
          
          calo_->geomInfo_.origin( CLHEP::Hep3Vector(xOrigin,0,0.5*(zCaloStart+zCaloEnd)) );          
          calo_->geomInfo_.trackerCenter(CLHEP::Hep3Vector(xOrigin,0,zTrackerCenter));
                    
          //Get the volume of the solid with hexagonal / rectangular base
          double hl              = calo_->caloInfo_.crystalZLength()/2.0;
          double ht              = calo_->caloInfo_.crystalXYLength()/2.0;
          double cryVolume       = 8*ht*ht*hl;
          calo_->caloInfo_.crystalVolume(cryVolume);

          //make sure above information has some degree of consistency, but run overlap checks to be sure
          CheckIt();

          // Create Disks
          MakeDisks();
    }


    DiskCalorimeterMaker::~DiskCalorimeterMaker() {}




    void DiskCalorimeterMaker::MakeDisks(void)
    {
        bool   crystalShift         = calo_->caloInfo_.crystalShift();
        double crystalHalfXY        = calo_->caloInfo_.crystalXYLength()/2.0;
        double crystalHalfZ         = calo_->caloInfo_.crystalZLength()/2.0;
        double wrapperHalfThick     = calo_->caloInfo_.wrapperThickness()/2.0;
        double wrapperHalfZ         = crystalHalfZ + wrapperHalfThick;
        double crystalCellRadius    = crystalHalfXY + 2.0*wrapperHalfThick;       
        double FPStyroHalfThick     = calo_->caloInfo_.FPStyrofoamZLength()/2.0;  
        double FPcarbonHalfThick    = calo_->caloInfo_.FPCarbonZLength()/2.0;  
        double FEEBoxHalfZ          = calo_->caloInfo_.FEEZLength()/2.0;
        double FEEBoxThick          = calo_->caloInfo_.FEEBoxThickness();
        double BPHoleHalfZ          = calo_->caloInfo_.BPHoleZLength()/2.0;       
        double innerCaseThickness   = calo_->caloInfo_.caseThicknessIn();
        double outerCaseThickness   = calo_->caloInfo_.caseThicknessOut();
        double outerRingEdgeThick   = calo_->caloInfo_.outerRingEdgeRLength();     
       
        double FPHalfZLength        = FPStyroHalfThick+2.0*FPcarbonHalfThick;
        double diskCaseHalfZLength  = wrapperHalfZ;
        double BPHalfZLength        = BPHoleHalfZ+FEEBoxHalfZ+2.0*FEEBoxThick;        

        double diskHalfZLength      = FPHalfZLength+diskCaseHalfZLength+BPHalfZLength;
        double motherHalfZ          = 0.5*(calo_->caloInfo_.envelopeZ1()-calo_->caloInfo_.envelopeZ0());         
        
        double crateFShieldThick    = calo_->caloInfo_.crateFShieldThickness();
        double FEBHalfZLength       = (calo_->caloInfo_.crateZLength()+calo_->caloInfo_.crateFShieldDeltaZ()+crateFShieldThick)/2.0;
                                              

        // Offsets between the disk and crystal cordinate systems (also disk and FF)
        CLHEP::Hep3Vector diskOriginToCrystalOrigin(0,0,-diskHalfZLength+2*FPHalfZLength+2.0*wrapperHalfThick);
        
        // Offset to align beginning of bottom shieSlding to beginning of crystal in z
        double crateToDiskDeltaZ    = FEBHalfZLength  -diskHalfZLength+2.0*FPHalfZLength;
       

        if (verbosityLevel_==99) std::cout<<"Disk components half length FP/Disk/BP "<<FPHalfZLength<<" "<<diskCaseHalfZLength<<" "<<BPHalfZLength<<std::endl;
        if (verbosityLevel_==99) std::cout<<"Disk diskOriginToCrystalOrigin "<<diskOriginToCrystalOrigin.z()<<std::endl;
        if (verbosityLevel_==99) std::cout<<"Disk crateToDiskDeltaZ "<<crateToDiskDeltaZ<<std::endl;
                
       
        for (int idisk=0; idisk<calo_->nDisks_; ++idisk)
        {
            double innerRadius = calo_->geomInfo_.diskCaseRadiusIn().at(idisk);
            double outerRadius = calo_->geomInfo_.diskCaseRadiusOut().at(idisk);
            double separation  = calo_->geomInfo_.diskZMotherShift().at(idisk);
           
            double dR1    = innerRadius - innerCaseThickness;
            double dR2    = outerRadius + outerCaseThickness + outerRingEdgeThick;
            double dZ     = 2.0*diskHalfZLength;

            CLHEP::Hep3Vector size(dR1,dR2,dZ) ;
            CLHEP::Hep3Vector originLocal(0, 0, -motherHalfZ + diskHalfZLength + separation);
            
	    CLHEP::Hep3Vector frontFaceCenter = calo_->geomInfo_.origin() + originLocal +diskOriginToCrystalOrigin;
            CLHEP::Hep3Vector backFaceCenter  = frontFaceCenter + CLHEP::Hep3Vector(0,0,2.0*diskCaseHalfZLength);


	    std::shared_ptr<Disk> thisDisk( new Disk(idisk,innerRadius, outerRadius, 2.0*crystalCellRadius, 
                                                     crystalShift, diskOriginToCrystalOrigin) );
            calo_->disks_.push_back(thisDisk);



            thisDisk->geomInfo().size(size);
            thisDisk->geomInfo().originLocal(originLocal);
            thisDisk->geomInfo().origin(calo_->geomInfo_.origin() + originLocal);
            thisDisk->geomInfo().rotation(CLHEP::HepRotation::IDENTITY*CLHEP::HepRotationZ( calo_->geomInfo().diskRotAngle().at(idisk)) );
            thisDisk->geomInfo().frontFaceCenter(frontFaceCenter);
            thisDisk->geomInfo().backFaceCenter(backFaceCenter);
            thisDisk->geomInfo().crateDeltaZ(crateToDiskDeltaZ);
	    thisDisk->geomInfo().envelopeRad(dR1,dR2);


            //fill the full Crystal List / diskId (direct access for performance optimization)
            int crystalOffset = calo_->fullCrystalList_.size();
            for (int icry=0;icry<thisDisk->nCrystals();++icry)
            {
                 Crystal& thisCrystal = thisDisk->crystal(icry);
                 calo_->fullCrystalList_.push_back(&thisCrystal);

                 //precompute the neighbors in the global frame
                 thisCrystal.setNeighbors(calo_->neighborsByLevel(icry+crystalOffset,1,false),false);
                 thisCrystal.setNeighbors(calo_->neighborsByLevel(icry+crystalOffset,1,true),true);
                 thisCrystal.setNextNeighbors(calo_->neighborsByLevel(icry+crystalOffset,2,false),false);
                 thisCrystal.setNextNeighbors(calo_->neighborsByLevel(icry+crystalOffset,2,true),true);

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
        double wrapperXYLength    = calo_->caloInfo_.crystalXYLength() + calo_->caloInfo_.wrapperThickness();      
        double wrapperZLength     = calo_->caloInfo_.crystalZLength() + calo_->caloInfo_.wrapperThickness();      
        double FPStyroZlength     = calo_->caloInfo_.FPStyrofoamZLength();  
        double FPCarbonZlength    = calo_->caloInfo_.FPCarbonZLength();  
        double FEEBoxZLength      = calo_->caloInfo_.FEEZLength();
        double BPHoleYLength      = calo_->caloInfo_.BPHoleYLength();
        double stripYLength       = calo_->caloInfo_.stripYLength();
        double BPHoleZLength      = calo_->caloInfo_.BPHoleZLength();
        double FEEBoxThick        = calo_->caloInfo_.FEEBoxThickness();        
        double innerCaseThickness = calo_->caloInfo_.caseThicknessIn();
        double outerCaseThickness = calo_->caloInfo_.caseThicknessOut();
        double outerRingEdgeThick = calo_->caloInfo_.outerRingEdgeRLength();     
       
        double FPZLength    = FPStyroZlength+2.0*FPCarbonZlength;
        double BPZLength    = FEEBoxZLength+BPHoleZLength+4.0*FEEBoxThick;
        double caseZLength  = wrapperZLength;

        double diskZLength  = FPZLength+caseZLength+BPZLength;
        

        //check number of readouts (now fixed to 2)
	if (nROPerCrystal !=2 && nROPerCrystal !=0) 
          throw cet::exception("DiskCaloGeom") << "calorimeter.readoutPerCrystal must be 2 at the moment (or 0 if no backplate is desired)\n";

        //check calorimeter fits inside mother envelope
        for (int i=0;i<calo_->nDisks_;++i)
        {
            if (calo_->geomInfo_.diskCaseRadiusIn().at(i) > calo_->geomInfo_.diskCaseRadiusOut().at(i))
                throw cet::exception("DiskCaloGeom") << "calorimeter.diskCaseRadiusIn > calorimeter.diskCaseRadiusOut for disk="<<i<<".\n";
            if (calo_->geomInfo_.diskCaseRadiusOut().at(i) + outerCaseThickness + outerRingEdgeThick > calo_->caloInfo_.envelopeRadiusOut())
                throw cet::exception("DiskCaloGeom") << "calorimeter outer radius larger than calorimeter mother for disk="<<i<<".\n";
            if (calo_->geomInfo_.diskCaseRadiusIn().at(i)- innerCaseThickness < calo_->caloInfo_.envelopeRadiusIn())
                throw cet::exception("DiskCaloGeom") << "calorimeter inner radius smaller than calorimeter mother for disk="<<i<<".\n";
            if (calo_->geomInfo_.diskZMotherShift().at(i) + diskZLength > calo_->caloInfo_.envelopeZ1())
                throw cet::exception("DiskCaloGeom") << "calorimeter length over mother envelope Z1 for disk="<<i<<".\n";
        }

        //check that holes in back plate are smaller than crystal, RO smaller than holes and FEE boxes fit
	if (nROPerCrystal)
        {
            if (calo_->caloInfo_.BPHoleXLength() > calo_->caloInfo_.crystalXYLength() || 
	         calo_->caloInfo_.BPHoleYLength() > calo_->caloInfo_.crystalXYLength() )
	         throw cet::exception("DiskCaloGeom") << "calorimeter backplate hole greater than crystal dimensions in X or Y \n";

	    if (calo_->caloInfo_.roXLength() > calo_->caloInfo_.BPHoleXLength()  || 
	         calo_->caloInfo_.roYLength() > calo_->caloInfo_.BPHoleYLength())
	         throw cet::exception("DiskCaloGeom") << "calorimeter readout larger than hole in X or Y \n";
	    if (calo_->caloInfo_.roZLength() > calo_->caloInfo_.BPHoleZLength())
	         throw cet::exception("DiskCaloGeom") << "calorimeter readout too thick to fit in hole \n";

	    if (calo_->caloInfo_.FEEXLength() > calo_->caloInfo_.BPHoleXLength()/nROPerCrystal)
	         throw cet::exception("DiskCaloGeom") << "calorimeter FEE box does not fit in X direction \n";
	    if (calo_->caloInfo_.FEEYLength() > calo_->caloInfo_.crystalXYLength()-2*FEEBoxThick)
	         throw cet::exception("DiskCaloGeom") << "calorimeter FEE box does not fit in Y direction \n";

            if (BPHoleYLength+stripYLength > wrapperXYLength )
       	         throw cet::exception("DiskCaloGeom") << "Box, strip and crystal length do not fit in Y direction \n";
        }
        
	//Check pipes
        for (int i=0;i<calo_->caloInfo().nPipes();++i)
        {
          if ( (calo_->caloInfo_.pipeTorRadius().at(i)- calo_->caloInfo_.pipeRadius()) <   calo_->caloInfo_.envelopeRadiusIn())
            throw cet::exception("DiskCaloGeom") << "element "<<i<<" of calorimeter.pipeTorRadius is smaller than disk inner radius\n";
          if ( (calo_->caloInfo_.pipeTorRadius().at(i)+ calo_->caloInfo_.pipeRadius()) >  calo_->caloInfo_.envelopeRadiusOut())
            throw cet::exception("DiskCaloGeom") << "element "<<i<<" of calorimeter.pipeTorRadius is larger than disk outer radius\n";
          if ( calo_->caloInfo_.pipeRadius()>  calo_->caloInfo_.FPStyrofoamZLength()/2.0)
            throw cet::exception("DiskCaloGeom") << "calorimeter pipe radius too large to fit inside styrofoam front panel\n";
        }
	
	
	//Just a few checks on crates
       int nBoards           = calo_->caloInfo_.nBoard();
       double radiatorDY     = calo_->caloInfo_.radiatorThickness()/2.0;
       double activeStripDY  = calo_->caloInfo_.activeStripThickness()/2.0;
       double passiveStripDY = calo_->caloInfo_.passiveStripThickness()/2.0;

       if ( nBoards*(radiatorDY+activeStripDY+passiveStripDY) > calo_->caloInfo_.crateYLength() ) 	
            throw cet::exception("DiskCaloGeom") << "calorimeter FEB boards too thick\n";	
       if (calo_->caloInfo_.crateFShieldYLength() > calo_->caloInfo_.crateYLength())
            throw cet::exception("DiskCaloGeom") << "calorimeter FEB front shile too long in Y direction\n";	    
       if (calo_->caloInfo_.crateSThickness() > calo_->caloInfo_.crateXLength())
            throw cet::exception("DiskCaloGeom") << "calorimeter FEB crate side too thick\n";
       if (calo_->caloInfo_.crateTThickness() > calo_->caloInfo_.crateYLength())
            throw cet::exception("DiskCaloGeom") << "calorimeter FEB crate top too thick\n";
       if (calo_->caloInfo_.crateXLength()/calo_->caloInfo_.crateRadiusIn() > calo_->caloInfo_.crateDeltaPhi()/180*3.141592)
            throw cet::exception("DiskCaloGeom") << "calorimeter FEB deltaPhi too small \n";
            

    }


}
