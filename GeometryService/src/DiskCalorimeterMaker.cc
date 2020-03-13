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


// For reference, git tag (ef94504f51edbbfeb54a5e63651856bdf5c0a60d) has generic placement of the disk origin.

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

#include <cmath>
#include <memory>



namespace mu2e {


    DiskCalorimeterMaker::DiskCalorimeterMaker(SimpleConfig const& config, double solenoidOffset)
    {
          
          calo_ = std::unique_ptr<DiskCalorimeter>(new DiskCalorimeter());
          
          //CALO INFO DATA
          verbosityLevel_ = config.getInt("calorimeter.verbosityLevel",0);

          calo_->nDisks_  = config.getInt("calorimeter.numberOfDisks");      

          calo_->caloInfo_.set("envelopeRadiusIn",       config.getDouble("calorimeter.caloMotherRadiusIn") );
          calo_->caloInfo_.set("envelopeRadiusOut",      config.getDouble("calorimeter.caloMotherRadiusOut") );
          calo_->caloInfo_.set("envelopeZ0",             config.getDouble("calorimeter.caloMotherZ0") );
          calo_->caloInfo_.set("envelopeZ1",             config.getDouble("calorimeter.caloMotherZ1") );
          calo_->caloInfo_.set("vdThickness",            config.getDouble("calorimeter.vdThickness") );
                    
          calo_->caloInfo_.set("diskInnerRingIn",        config.getDouble("calorimeter.diskInnerRingIn") );
          calo_->caloInfo_.set("diskInnerRingOut",       config.getDouble("calorimeter.diskInnerRingOut") );
          calo_->caloInfo_.set("diskCrystalRadiusIn",    config.getDouble("calorimeter.diskCrystalRadiusIn") );
          calo_->caloInfo_.set("diskCrystalRadiusOut",   config.getDouble("calorimeter.diskCrystalRadiusOut") );          
          calo_->caloInfo_.set("diskOuterRingIn",        config.getDouble("calorimeter.diskOuterRingIn") );
          calo_->caloInfo_.set("diskOuterRingOut",       config.getDouble("calorimeter.diskOuterRingOut") );
          calo_->caloInfo_.set("diskCaseZLength",        config.getDouble("calorimeter.diskCaseZLength") );
	  calo_->caloInfo_.set("diskOutRingEdgeZLength", config.getDouble("calorimeter.diskOutRingEdgeZLength") );
	  calo_->caloInfo_.set("diskOutRingEdgeRLength", config.getDouble("calorimeter.diskOutRingEdgeRLength") );
	  calo_->caloInfo_.set("diskStepThickness",      config.getDouble("calorimeter.diskStepThickness") );
          std::vector<double> temp;
          config.getVectorDouble("calorimeter.diskZ0MotherShift", temp, calo_->nDisks_);
          calo_->caloInfo_.set("diskZMotherShift",temp);

	  calo_->caloInfo_.set("crystalXYLength",        config.getDouble("calorimeter.crystalXYLength") );
	  calo_->caloInfo_.set("crystalZLength",         config.getDouble("calorimeter.crystalZLength") );
	  calo_->caloInfo_.set("crystalFrameZLength",    config.getDouble("calorimeter.crystalFrameZLength") );
	  calo_->caloInfo_.set("crystalFrameThickness",  config.getDouble("calorimeter.crystalFrameThickness") );
	  calo_->caloInfo_.set("wrapperThickness",       config.getDouble("calorimeter.wrapperThickness") );
	  calo_->caloInfo_.set("refractiveIndex",        config.getDouble("calorimeter.refractiveIndex") );	  
	  calo_->caloInfo_.set("readoutPerCrystal",      config.getInt("calorimeter.readoutPerCrystal") );
	  calo_->caloInfo_.set("readoutXLength",         config.getDouble("calorimeter.readoutXLength") );
	  calo_->caloInfo_.set("readoutYLength",         config.getDouble("calorimeter.readoutYLength") );
	  calo_->caloInfo_.set("readoutZLength",         config.getDouble("calorimeter.readoutZLength") );
	            
          calo_->caloInfo_.set("FEEXLength",             config.getDouble("calorimeter.FEEXLength") );
	  calo_->caloInfo_.set("FEEYLength",             config.getDouble("calorimeter.FEEYLength") );
	  calo_->caloInfo_.set("FEEZLength",             config.getDouble("calorimeter.FEEZLength") );
	  calo_->caloInfo_.set("FEEBoxThickness",        config.getDouble("calorimeter.FEEBoxThickness") );
	  calo_->caloInfo_.set("BPStripThickness",       config.getDouble("calorimeter.BPStripThickness") );
	  calo_->caloInfo_.set("BPHoleXLength",          config.getDouble("calorimeter.BPHoleXLength") );
	  calo_->caloInfo_.set("BPHoleYLength",          config.getDouble("calorimeter.BPHoleYLength") );
	  calo_->caloInfo_.set("BPHoleZLength",          config.getDouble("calorimeter.BPHoleZLength") );
	  calo_->caloInfo_.set("BPOuterRadius",          config.getDouble("calorimeter.BPOuterRadius") );
	  calo_->caloInfo_.set("BPPipeRadiusHigh",       config.getDouble("calorimeter.BPPipeRadiusHigh") );
	  calo_->caloInfo_.set("BPPipeRadiusLow",        config.getDouble("calorimeter.BPPipeRadiusLow") );
	  calo_->caloInfo_.set("BPPipeThickness",        config.getDouble("calorimeter.BPPipeThickness") );
	  calo_->caloInfo_.set("BPPipeZOffset",          config.getDouble("calorimeter.BPPipeZOffset") );
	  
	  calo_->caloInfo_.set("FPInnerRadius",          config.getDouble("calorimeter.FPInnerRadius") );
	  calo_->caloInfo_.set("FPOuterRadius",          config.getDouble("calorimeter.FPOuterRadius") );
	  calo_->caloInfo_.set("FPFoamZLength",          config.getDouble("calorimeter.FPFoamZLength") );
	  calo_->caloInfo_.set("FPCarbonZLength",        config.getDouble("calorimeter.FPCarbonZLength") );	  
          calo_->caloInfo_.set("FPCoolPipeTorRadius",    config.getDouble("calorimeter.FPCoolPipeTorRadius") );
          calo_->caloInfo_.set("FPCoolPipeRadius",       config.getDouble("calorimeter.FPCoolPipeRadius") );
          calo_->caloInfo_.set("FPCoolPipeThickness",    config.getDouble("calorimeter.FPCoolPipeThickness") );
	  calo_->caloInfo_.set("nPipes",                 config.getInt("calorimeter.nPipes") );
	  calo_->caloInfo_.set("pipeRadius",             config.getDouble("calorimeter.pipeRadius") );
	  calo_->caloInfo_.set("pipeThickness",          config.getDouble("calorimeter.pipeThickness") );
	  calo_->caloInfo_.set("pipeInitSeparation",     config.getDouble("calorimeter.pipeInitSeparation") );
          temp.clear();
            config.getVectorDouble("calorimeter.pipeTorRadius", temp, calo_->caloInfo_.getInt("nPipes"));
	  calo_->caloInfo_.set("pipeTorRadius",temp );

	  calo_->caloInfo_.set("numberOfCrates",         config.getInt("calorimeter.numberOfCrates") );
	  calo_->caloInfo_.set("nCrateBeforeSpace",      config.getInt("calorimeter.nCrateBeforeSpace") );
	  calo_->caloInfo_.set("numberOfBoards",         config.getInt("calorimeter.numberOfBoards") );
	  calo_->caloInfo_.set("crateXLength",           config.getDouble("calorimeter.crateXLength") );
	  calo_->caloInfo_.set("crateYLength",           config.getDouble("calorimeter.crateYLength") );
	  calo_->caloInfo_.set("crateZLength",           config.getDouble("calorimeter.crateZLength") );
	  calo_->caloInfo_.set("crateFShieldThickness",  config.getDouble("calorimeter.crateFShieldThickness") );
	  calo_->caloInfo_.set("crateBShieldThickness",  config.getDouble("calorimeter.crateBShieldThickness") );
	  calo_->caloInfo_.set("crateTThickness",        config.getDouble("calorimeter.crateTThickness") );
	  calo_->caloInfo_.set("crateSThickness",        config.getDouble("calorimeter.crateSThickness") );
	  calo_->caloInfo_.set("crateFShieldYLength",    config.getDouble("calorimeter.crateFShieldYLength") );
	  calo_->caloInfo_.set("crateFShieldDeltaZ",     config.getDouble("calorimeter.crateFShieldDeltaZ") );
	  calo_->caloInfo_.set("cratephi0",              config.getDouble("calorimeter.cratephi0") );
	  calo_->caloInfo_.set("crateDeltaPhi",          config.getDouble("calorimeter.crateDeltaPhi") );
	  calo_->caloInfo_.set("radiatorThickness",      config.getDouble("calorimeter.radiatorThickness") );
	  calo_->caloInfo_.set("radiatorZLength",        config.getDouble("calorimeter.radiatorZLength") );
	  calo_->caloInfo_.set("activeStripThickness",   config.getDouble("calorimeter.activeStripThickness") );
	  calo_->caloInfo_.set("passiveStripThickness",  config.getDouble("calorimeter.passiveStripThickness") );
	  	 	  
          
          if (calo_->caloInfo_.getInt("nPipes")==0)
          {
             calo_->caloInfo_.set("FPInnerRadius",0.0);
             calo_->caloInfo_.set("FPOuterRadius",0.0);
             calo_->caloInfo_.set("FPFoamZLength",0.0);
             calo_->caloInfo_.set("FPCarbonZLength",0.0);
             calo_->caloInfo_.set("FPCoolPipeRadius",0.0);
             calo_->caloInfo_.set("FPCoolPipeTorRadius",0.0);
             calo_->caloInfo_.set("pipeRadius",0.0);         
             calo_->caloInfo_.set("pipeThickness",0.0);         
             calo_->caloInfo_.set("pipeInitSeparation",0.0);         
          }
          
          if (calo_->caloInfo_.getInt("readoutPerCrystal")==0)
          {
	      calo_->caloInfo_.set("FEEZLength",0.0 );
 	      calo_->caloInfo_.set("BPStripThickness",0.0 );
	      calo_->caloInfo_.set("FEEBoxThickness",0.0 );
	      calo_->caloInfo_.set("BPHoleZLength",0.0 );
 	      calo_->caloInfo_.set("readoutZLength",0.0 );
          }

          
          //CALORIMETER ORIGIN AND FRONT FACE (FF)
          double zTrackerCenter  =  config.getDouble("mu2e.detectorSystemZ0");
          double xOrigin         = -config.getDouble("mu2e.solenoidOffset");
          double zCaloStart      =  config.getDouble("calorimeter.caloMotherZ0");
          double zCaloEnd        =  config.getDouble("calorimeter.caloMotherZ1");
          
          calo_->geomUtil_.origin( CLHEP::Hep3Vector(xOrigin,0,0.5*(zCaloStart+zCaloEnd)) );          
          calo_->geomUtil_.trackerCenter(CLHEP::Hep3Vector(xOrigin,0,zTrackerCenter));
          calo_->geomUtil_.crystalZLength(config.getDouble("calorimeter.crystalZLength"));


	  // CACHE THIS ONE FOR EFFICIENCY (REALLY NEEDED SO DON'T REMOVE)
          calo_->caloInfo_.nROPerCrystal( config.getInt(   "calorimeter.readoutPerCrystal") );



          // CALCULATE THE TOTAL Z LENGTH OF THE FULL DISK AND SUB-COMPONENTS 
          // LOOK AT CONSTRUCTDISKCALORIMETER TO GET THE APPROPRIATE FORMULAS FOR 
          // frontPanelHalfThick, diskCaseHalfZLength, BPFEEHalfZ and crateFullDZ 
          //
          double vdThickness             = calo_->caloInfo_.getDouble("vdThickness");  
          double FPCarbonThick           = calo_->caloInfo_.getDouble("FPCarbonZLength");  
          double FPFoamThick             = calo_->caloInfo_.getDouble("FPFoamZLength");  
          double FPCoolPipeRadius        = calo_->caloInfo_.getDouble("FPCoolPipeRadius");  
          double FPpipeRadius            = calo_->caloInfo_.getDouble("pipeRadius");
          double diskCaseHalfZLength     = calo_->caloInfo_.getDouble("diskCaseZLength")/2.0;
          double crystalHalfZLength      = calo_->caloInfo_.getDouble("crystalZLength")/2.0;
          double crystalFrameHalfZLength = calo_->caloInfo_.getDouble("crystalFrameZLength")/2.0;
          double BPHoleHalfZ             = calo_->caloInfo_.getDouble("BPHoleZLength")/2.0;       
          double FEEBoxHalfZ             = calo_->caloInfo_.getDouble("FEEZLength")/2.0;
          double FEEBoxThick             = calo_->caloInfo_.getDouble("FEEBoxThickness");
          double BPPipeRadiusHigh        = calo_->caloInfo_.getDouble("BPPipeRadiusHigh");
          double BPPipeHalfZOffset       = calo_->caloInfo_.getDouble("BPPipeZOffset")/2.0;
          double crateZLength            = calo_->caloInfo_.getDouble("crateZLength");
          double crateFShieldThick       = calo_->caloInfo_.getDouble("crateFShieldThickness");
          double crateFShieldDeltaZ      = calo_->caloInfo_.getDouble("crateFShieldDeltaZ");
                   
          FPHalfZLength_        = (FPCarbonThick+FPFoamThick-FPpipeRadius+FPCoolPipeRadius)/2.0;
          diskCaseHalfZLength_  = diskCaseHalfZLength;        
          BPHalfZLength_        = BPHoleHalfZ+FEEBoxHalfZ+2.0*FEEBoxThick+BPPipeHalfZOffset+BPPipeRadiusHigh;        
          diskHalfZLength_      = FPHalfZLength_+diskCaseHalfZLength_+BPHalfZLength_ + vdThickness;
          FEBHalfZLength_       = (crateZLength+crateFShieldDeltaZ+crateFShieldThick + 2.0*vdThickness)/2.0;
          motherHalfZ_          = (calo_->caloInfo_.getDouble("envelopeZ1")-calo_->caloInfo_.getDouble("envelopeZ0"))/2.0;

          // OFFSET TO ALIGN BEGINNING OF BOTTOM SHIESLDING TO BEGINNING OF CRYSTAL IN Z
          crateToDiskDeltaZ_ = FEBHalfZLength_ - diskHalfZLength_ + 2.0*FPHalfZLength_ + vdThickness;

    

      
          // OFFSETS BETWEEN THE DISK AND CRYSTAL CORDINATE SYSTEMS, I.E. DISTANCE BETWEEN FRONT FACE DISK AND FRONT FACE CRYSTALS
          // LOOK AT CRYSTALPOSITION
          double disp = -diskHalfZLength_+vdThickness + 2*FPHalfZLength_+ 2.0*(diskCaseHalfZLength - crystalHalfZLength - crystalFrameHalfZLength);
          diskOriginToCrystalOrigin_ = CLHEP::Hep3Vector(0,0,disp);
          
         if (verbosityLevel_==99) std::cout<<"Disk components half length FP/Disk/BP "<<FPHalfZLength_<<" "<<diskCaseHalfZLength_<<" "<<BPHalfZLength_<<std::endl;
         if (verbosityLevel_==99) std::cout<<"Disk diskOriginToCrystalOrigin_ "<<diskOriginToCrystalOrigin_.z()<<std::endl;
         if (verbosityLevel_==99) std::cout<<"Disk crateToDiskDeltaZ_ "<<crateToDiskDeltaZ_<<std::endl;
                    
         CheckIt();          
         MakeIt();
    }


    DiskCalorimeterMaker::~DiskCalorimeterMaker() {}



    void DiskCalorimeterMaker::MakeIt(void)
    {
        
        double crystalHalfXY        = calo_->caloInfo_.getDouble("crystalXYLength")/2.0;
        double wrapperHalfThick     = calo_->caloInfo_.getDouble("wrapperThickness")/2.0;
        double crystalCellRadius    = crystalHalfXY + 2.0*wrapperHalfThick;       
        double innerCaseRadius      = calo_->caloInfo_.getDouble("diskInnerRingIn");
        double outerCaseRadius      = calo_->caloInfo_.getDouble("diskOuterRingOut");
        double innerCrysRadius      = calo_->caloInfo_.getDouble("diskCrystalRadiusIn");
        double outerCrysRadius      = calo_->caloInfo_.getDouble("diskCrystalRadiusOut");
        double outerRingEdgeThick   = calo_->caloInfo_.getDouble("diskOutRingEdgeRLength");                                       
                      
       
        for (int idisk=0; idisk<calo_->nDisks_; ++idisk)
        {
            int crystalOffset      = calo_->fullCrystalList_.size();
            
            double separation      = calo_->caloInfo_.getVDouble("diskZMotherShift").at(idisk);
            double angleZ          = 0;
           
            double dR1 = innerCaseRadius;
            double dR2 = outerCaseRadius + outerRingEdgeThick;
            double dZ  = 2.0*diskHalfZLength_;
            
            CLHEP::Hep3Vector size(dR1,dR2,dZ) ;
            CLHEP::Hep3Vector originLocal(0, 0, -motherHalfZ_ + diskHalfZLength_ + separation);
            
	    CLHEP::Hep3Vector frontFaceCenter = calo_->geomUtil_.origin() + originLocal +diskOriginToCrystalOrigin_;
            CLHEP::Hep3Vector backFaceCenter  = frontFaceCenter + CLHEP::Hep3Vector(0,0,2.0*diskCaseHalfZLength_);
            CLHEP::HepRotation diskRotation   = CLHEP::HepRotation::IDENTITY*CLHEP::HepRotationZ(angleZ);

	    std::shared_ptr<Disk> thisDisk( new Disk(idisk,innerCaseRadius, outerCaseRadius,innerCrysRadius,outerCrysRadius, 
                                                     2.0*crystalCellRadius, crystalOffset, diskOriginToCrystalOrigin_) );
            calo_->disks_.push_back(thisDisk);

            thisDisk->geomInfo().size(size);
            thisDisk->geomInfo().originLocal(originLocal);
            thisDisk->geomInfo().origin(calo_->geomUtil_.origin() + originLocal);
            thisDisk->geomInfo().rotation(diskRotation);
            thisDisk->geomInfo().frontFaceCenter(frontFaceCenter);
            thisDisk->geomInfo().backFaceCenter(backFaceCenter);
            thisDisk->geomInfo().crateDeltaZ(crateToDiskDeltaZ_);
	    thisDisk->geomInfo().envelopeRad(dR1,dR2);


            //fill the full Crystal List / diskId (direct access for performance optimization)
            for (unsigned icry=0;icry<thisDisk->nCrystals();++icry)
            {
                 Crystal& thisCrystal = thisDisk->crystal(icry);
                 calo_->fullCrystalList_.push_back(&thisCrystal);

                 //precompute the neighbors in the global frame
                 thisCrystal.setNeighbors(calo_->neighborsByLevel(icry + crystalOffset,1,false),false);
                 thisCrystal.setNeighbors(calo_->neighborsByLevel(icry + crystalOffset,1,true),true);
                 thisCrystal.setNextNeighbors(calo_->neighborsByLevel(icry + crystalOffset,2,false),false);
                 thisCrystal.setNextNeighbors(calo_->neighborsByLevel(icry + crystalOffset,2,true),true);

                 //pre-compute the crystal position in the mu2e frame (aka global frame)
                 CLHEP::Hep3Vector globalPosition = thisDisk->geomInfo().origin() + thisDisk->geomInfo().inverseRotation()*(thisCrystal.localPosition());
                 thisCrystal.setPosition(globalPosition);                                  
            }
            
            //calculate the position of the inner and outer steps. 
            std::vector<double> stepsInX, stepsInY, stepsOutX, stepsOutY;
            
            int nrowIn = int(innerCaseRadius/2.0/crystalCellRadius)+1;
            for (int irow=-nrowIn;irow<=nrowIn;++irow)
            {  
                int icry = thisDisk->idMinCrystalInside(irow);
                if (icry <0 || icry > int(thisDisk->nCrystals())) continue;
                double xmin = thisDisk->crystal(icry).localPosition().x()-crystalCellRadius;
                double ymin = thisDisk->crystal(icry).localPosition().y();
                if (xmin > crystalCellRadius) {stepsInX.push_back(xmin);stepsInY.push_back(ymin);}               
            }

            int nrowOut = int(outerCaseRadius/2.0/crystalCellRadius)+1;
            for (int irow=-nrowOut;irow<=nrowOut;++irow)
            {  
                int icry = thisDisk->idMaxCrystalInside(irow);
                if (icry <0 || icry > int(thisDisk->nCrystals())) continue;
                double xmax = thisDisk->crystal(icry).localPosition().x()+crystalCellRadius;
                double ymax = thisDisk->crystal(icry).localPosition().y();
                if (xmax > crystalCellRadius) {stepsOutX.push_back(xmax);stepsOutY.push_back(ymax);}               
            }
            
            calo_->caloInfo_.set("stepsInsideX",stepsInX);
            calo_->caloInfo_.set("stepsInsideY",stepsInY);            
            calo_->caloInfo_.set("stepsOutsideX",stepsOutX);
            calo_->caloInfo_.set("stepsOutsideY",stepsOutY);

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
        //check number of readouts (now fixed to 2)
        int nROPerCrystal = calo_->caloInfo_.getInt("readoutPerCrystal");
	if (nROPerCrystal !=2 && nROPerCrystal !=0) 
          throw cet::exception("DiskCaloGeom") << "calorimeter.readoutPerCrystal must be 2 at the moment (or 0 if no backplate is desired)\n";


        //check calorimeter fits inside mother envelope
        double diskRin            = calo_->caloInfo_.getDouble("diskInnerRingIn");
        double diskRout           = calo_->caloInfo_.getDouble("diskOuterRingOut");
        double outerRingEdgeThick = calo_->caloInfo_.getDouble("diskOutRingEdgeRLength");     
        for (int i=0;i<calo_->nDisks_;++i)
        {            
            if (diskRin > diskRout)
                  throw cet::exception("DiskCaloGeom") << "calorimeter.diskInnerRingIn > calorimeter.diskOuterRingOut for disk="<<i<<".\n";
            
            if (diskRout + outerRingEdgeThick > calo_->caloInfo_.getDouble("envelopeRadiusOut"))
                  throw cet::exception("DiskCaloGeom") << "calorimeter outer radius larger than calorimeter mother for disk="<<i<<".\n";
            
            if (diskRin < calo_->caloInfo_.getDouble("envelopeRadiusIn"))
                  throw cet::exception("DiskCaloGeom") << "calorimeter inner radius smaller than calorimeter mother for disk="<<i<<".\n";
            
            if (calo_->caloInfo_.getVDouble("diskZMotherShift").at(i) + 2*diskHalfZLength_ > calo_->caloInfo_.getDouble("envelopeZ1"))
                  throw cet::exception("DiskCaloGeom") << "calorimeter length over mother envelope Z1 for disk="<<i<<".\n";
        }

        //check that holes in back plate are smaller than crystal, RO smaller than holes and FEE boxes fit
	if (nROPerCrystal)
        {
            if (calo_->caloInfo_.getDouble("BPHoleXLength") > calo_->caloInfo_.getDouble("crystalXYLength") || 
	        calo_->caloInfo_.getDouble("BPHoleYLength") > calo_->caloInfo_.getDouble("crystalXYLength") )
	          throw cet::exception("DiskCaloGeom") << "calorimeter backplate hole greater than crystal dimensions in X or Y \n";

	    if (calo_->caloInfo_.getDouble("readoutXLength") > calo_->caloInfo_.getDouble("BPHoleXLength")  || 
	        calo_->caloInfo_.getDouble("readoutYLength") > calo_->caloInfo_.getDouble("BPHoleYLength"))
	          throw cet::exception("DiskCaloGeom") << "calorimeter readout larger than hole in X or Y \n";
	    
            if (calo_->caloInfo_.getDouble("readoutZLength") > calo_->caloInfo_.getDouble("BPHoleZLength"))
	          throw cet::exception("DiskCaloGeom") << "calorimeter readout too thick to fit in hole \n"; 

	    if (calo_->caloInfo_.getDouble("FEEXLength") > calo_->caloInfo_.getDouble("BPHoleXLength")/nROPerCrystal)
	          throw cet::exception("DiskCaloGeom") << "calorimeter FEE box does not fit in X direction \n";
	    
            if (calo_->caloInfo_.getDouble("FEEYLength") > calo_->caloInfo_.getDouble("crystalXYLength")-2*calo_->caloInfo_.getDouble("FEEBoxThickness"))
	          throw cet::exception("DiskCaloGeom") << "calorimeter FEE box does not fit in Y direction \n";
        }

	//Check pipes       
        for (int i=0;i<calo_->caloInfo().getInt("nPipes");++i)
        {
           double pipeTorRadius     = calo_->caloInfo_.getVDouble("pipeTorRadius").at(i);
           double pipeRadius        = calo_->caloInfo_.getDouble("pipeRadius");
           double radiusIn          = calo_->caloInfo_.getDouble("FPInnerRadius");
           double radiusOut         = calo_->caloInfo_.getDouble("FPOuterRadius");
           double foamZLength       = calo_->caloInfo_.getDouble("FPFoamZLength");
           double carbonThick       = calo_->caloInfo_.getDouble("FPCarbonZLength");  
           double foamThick         = calo_->caloInfo_.getDouble("FPFoamZLength");  
           double coolPipeTorRadius = calo_->caloInfo_.getDouble("FPCoolPipeTorRadius");  
           double coolPipeRadius    = calo_->caloInfo_.getDouble("FPCoolPipeRadius");  
           
           if ( pipeTorRadius - pipeRadius < radiusIn)
                 throw cet::exception("DiskCaloGeom") << "element "<<i<<" of calorimeter.pipeTorRadius is smaller than disk inner radius\n";
           
           if ( pipeTorRadius + pipeRadius > radiusOut)
                 throw cet::exception("DiskCaloGeom") << "element "<<i<<" of calorimeter.pipeTorRadius is larger than disk outer radius\n";
           
           if ( pipeRadius > foamZLength/2.0)
                 throw cet::exception("DiskCaloGeom") << "calorimeter pipe radius too large to fit inside Foam front panel\n";
                      
           if ( carbonThick + foamThick - pipeRadius < coolPipeRadius)
                 throw cet::exception("DiskCaloGeom") << "calorimeter cooling pipe radius too large\n";
           
           if ( carbonThick + pipeRadius > coolPipeRadius)
                 throw cet::exception("DiskCaloGeom") << "calorimeter cooling pipe radius too small\n";
           
           if ( coolPipeTorRadius - coolPipeRadius < radiusOut)
                 throw cet::exception("DiskCaloGeom") << "cooling pipe too large, overlap with foam structure\n";                 
        }
	
	
	//Just a few checks on crates
        int nBoards           = calo_->caloInfo_.getInt("numberOfBoards");
        double radiatorDY     = calo_->caloInfo_.getDouble("radiatorThickness")/2.0;
        double activeStripDY  = calo_->caloInfo_.getDouble("activeStripThickness")/2.0;
        double passiveStripDY = calo_->caloInfo_.getDouble("passiveStripThickness")/2.0;
        double crateXLength   = calo_->caloInfo_.getDouble("crateXLength");
        double crateYLength   = calo_->caloInfo_.getDouble("crateYLength");

        if ( nBoards*(radiatorDY+activeStripDY+passiveStripDY) > crateYLength) 	
              throw cet::exception("DiskCaloGeom") << "calorimeter FEB boards too thick\n";	
        
        if (calo_->caloInfo_.getDouble("crateFShieldYLength") > crateYLength)
              throw cet::exception("DiskCaloGeom") << "calorimeter FEB front shile too long in Y direction\n";	    
        
        if (calo_->caloInfo_.getDouble("crateSThickness") > crateXLength)
              throw cet::exception("DiskCaloGeom") << "calorimeter FEB crate side too thick\n";
        
        if (calo_->caloInfo_.getDouble("crateTThickness") > crateYLength )
              throw cet::exception("DiskCaloGeom") << "calorimeter FEB crate top too thick\n";            

    }


}
