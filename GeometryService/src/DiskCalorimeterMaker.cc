//
// Make a Disk Calorimeter.
//
// original authors  Bertrand Echenarrd
//
// Disk geometry
//
//  see Mu2eG4/src/constructDiskCalorimeter.cc for the geometry
//    the disk are split in three pieces: front palte, disk case and back plate
//
//  front plate has calibration pipes
//  middle case has crystal enveloped in wrapper with front cap (back of the crystal is free)
//  back plate has a solid plate with holes and readouts inside + FEE box outside
//
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
//

#include "cetlib_except/exception.h"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CalorimeterGeom/inc/Disk.hh"
#include "Offline/CalorimeterGeom/inc/Crystal.hh"
#include "Offline/GeometryService/inc/DiskCalorimeterMaker.hh"

#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Vector/ThreeVector.h"




namespace mu2e {

    DiskCalorimeterMaker::DiskCalorimeterMaker(SimpleConfig const& config, double solenoidOffset)
    {

       calo_ = std::unique_ptr<DiskCalorimeter>(new DiskCalorimeter());

       verbosityLevel_ = config.getInt("calorimeter.verbosityLevel",0);

       calo_->G4Info_.set("caloDiskRadiusIn",       config.getDouble("calorimeter.caloDiskRadiusIn") );    //inner radius of disk enveloppe
       calo_->G4Info_.set("caloDiskRadiusOut",      config.getDouble("calorimeter.caloDiskRadiusOut") );   //outer radius of disk enveloppe
       calo_->G4Info_.set("caloFEBRadiusOut",       config.getDouble("calorimeter.caloFEBRadiusOut") );    //outer radius of FEB enveloppe
       calo_->G4Info_.set("caloMotherZ0",           config.getDouble("calorimeter.caloMotherZ0") );        //upstream Z of calorimeter enveloppe volume
       calo_->G4Info_.set("caloMotherZ1",           config.getDouble("calorimeter.caloMotherZ1") );        //downstream Z of calorimeter enveloppe volume
       calo_->G4Info_.set("caloMotherXorig",        config.getDouble("mu2e.solenoidOffset") );             //x coordinate center of calo mother volume
       calo_->G4Info_.set("vdThickness",            config.getDouble("calorimeter.vdThickness") );         //virtual detector thickness

       calo_->G4Info_.set("diskInAlRingRIn",        config.getDouble("calorimeter.diskInAlRingRIn") );     //inner radius of Al ring inside disk
       calo_->G4Info_.set("diskInAlRingZLength",    config.getDouble("calorimeter.diskInAlRingZLength") ); //Z length of Al ring inside disk
       calo_->G4Info_.set("diskInCFRingRIn",        config.getDouble("calorimeter.diskInCFRingRIn") );     //inner radius of carbon fiber ring
       calo_->G4Info_.set("diskInCFRingROut",       config.getDouble("calorimeter.diskInCFRingROut") );    //outer radius of carbon fiber ring
       calo_->G4Info_.set("diskCaseZLength",        config.getDouble("calorimeter.crystalZLength")         //Z Length of crystal support structure  (large Al ring)
                                                    + config.getDouble("calorimeter.crystalCapZLength") );
       calo_->G4Info_.set("diskCrystalRIn",         config.getDouble("calorimeter.diskCrystalRIn") );      //inner radius of volume containing crystals
       calo_->G4Info_.set("diskCrystalROut",        config.getDouble("calorimeter.diskCrystalROut") );     //outer radius of volume containing crystals
       calo_->G4Info_.set("diskCaseRingROut",       config.getDouble("calorimeter.diskCaseRingROut") );    //outer radius of crystal support structure (large Al ring)
       calo_->G4Info_.set("diskOutRailZLength",     config.getDouble("calorimeter.diskOutRailZLength") );  //Z length of rail outside the crystal support structure
       calo_->G4Info_.set("diskOutRailROut",        config.getDouble("calorimeter.diskOutRailROut") );     //outer radius of rail outside the crystal support structure
       calo_->G4Info_.set("diskStepThickness",      config.getDouble("calorimeter.diskStepThickness") );   //thickness of carbon fiber steps glued to inner hole
       calo_->G4Info_.set("diskZMotherShift",       getVDouble(config,"diskZMotherShift",CaloConst::_nDisk)); //Z offset betweenfront face of disk+FEB envelope and calo mother volume Z0
       calo_->G4Info_.set("diskCrystalFile",        config.getString("calorimeter.diskCrystalFile"));      //File with measured location of crystals

       calo_->G4Info_.set("crystalXYLength",        config.getDouble("calorimeter.crystalXYLength") );     //crysal transverse dimension
       calo_->G4Info_.set("crystalZLength",         config.getDouble("calorimeter.crystalZLength") );      //crystal Z length
       calo_->G4Info_.set("crystalCapZLength",      config.getDouble("calorimeter.crystalCapZLength") );   //crystal cap Z length (only in front of crystal)
       calo_->G4Info_.set("wrapperThickness",       config.getDouble("calorimeter.wrapperThickness") );    //crystal wrapper thickness
       calo_->G4Info_.set("refractiveIndex",        config.getDouble("calorimeter.refractiveIndex") );     //crystal refractive index

       calo_->G4Info_.set("nSiPMPerCrystal",        int(CaloConst::_nSiPMPerCrystal));                     //number of SiPM per crystal
       calo_->G4Info_.set("readoutXLength",         config.getDouble("calorimeter.readoutXLength") );      //SiPM x length
       calo_->G4Info_.set("readoutYLength",         config.getDouble("calorimeter.readoutYLength") );      //SiPM Y length
       calo_->G4Info_.set("readoutZLength",         config.getDouble("calorimeter.readoutZLength") );      //SiPM Z length
       calo_->G4Info_.set("shimStepsInRowId",       getVInt(config,"shimStepsInRowId") );                  //indices of crystal rows w.r.t center calorimeter with inner shims
       calo_->G4Info_.set("shimStepsOutRowId",      getVInt(config,"shimStepsOutRowId") );                 //indices of crystal rows w.r.t center calorimeter with outer shims
       std::vector<int> tempInt;
       for (uint16_t i=0;i<CaloConst::_nCaphriCrystal;++i) tempInt.push_back(CaloConst::_caphriId[i]);
       calo_->G4Info_.set("caphriCrystalId",tempInt);                                                      //LYSO crystal ids

       calo_->G4Info_.set("hasCrates",              config.getBool  ("calorimeter.hasCrates") );           //include or exclude crates
       calo_->G4Info_.set("hasFrontPanel",          config.getBool  ("calorimeter.hasFrontPanel") );       //include or exclude front panel with calibration pipes
       calo_->G4Info_.set("hasBackPanel",           config.getBool  ("calorimeter.hasBackPanel") );        //include or exclude back panel with readouts and cooling pipes
       calo_->G4Info_.set("FEEXLength",             config.getDouble("calorimeter.FEEXLength") );          //X length of FEE copper box
       calo_->G4Info_.set("FEEYLength",             config.getDouble("calorimeter.FEEYLength") );          //Y length of FEE copper box
       calo_->G4Info_.set("FEEZLength",             config.getDouble("calorimeter.FEEZLength") );          //Z length of FEE copper box
       calo_->G4Info_.set("FEEBoxThickness",        config.getDouble("calorimeter.FEEBoxThickness") );     //FEE copper box thickness
       calo_->G4Info_.set("BPStripThickness",       config.getDouble("calorimeter.BPStripThickness") );    //Back plate cooling strip thickness
       calo_->G4Info_.set("BPHoleXLength",          config.getDouble("calorimeter.BPHoleXLength") );       //X length of hole in back plate
       calo_->G4Info_.set("BPHoleYLength",          config.getDouble("calorimeter.BPHoleYLength") );       //Y length of hole in back plate
       calo_->G4Info_.set("BPHoleZLength",          config.getDouble("calorimeter.BPHoleZLength") );       //Z length of hole in back plate
       calo_->G4Info_.set("BPOuterRadius",          config.getDouble("calorimeter.BPOuterRadius") );       //Back plate outer radius
       calo_->G4Info_.set("BPPipeRadiusHigh",       config.getDouble("calorimeter.BPPipeRadiusHigh") );    //Radius of large back plate outer pipe
       calo_->G4Info_.set("BPPipeRadiusLow",        config.getDouble("calorimeter.BPPipeRadiusLow") );     //Radius of small back plate outer pipe
       calo_->G4Info_.set("BPPipeThickness",        config.getDouble("calorimeter.BPPipeThickness") );     //Thickness of back plate outer pipe
       calo_->G4Info_.set("BPPipeZOffset",          config.getDouble("calorimeter.BPPipeZOffset") );       //Z offset between outer pipe center and back plate

       calo_->G4Info_.set("FPInnerRadius",          config.getDouble("calorimeter.FPInnerRadius") );       //Front plate inner radius
       calo_->G4Info_.set("FPOuterRadius",          config.getDouble("calorimeter.FPOuterRadius") );       //Front plate outer radius
       calo_->G4Info_.set("FPFoamZLength",          config.getDouble("calorimeter.FPFoamZLength") );       //Thickness of foam inside front plate
       calo_->G4Info_.set("FPCarbonZLength",        config.getDouble("calorimeter.FPCarbonZLength") );     //Z length of carbon fiber panel on front plate
       calo_->G4Info_.set("FPCoolPipeTorRadius",    config.getDouble("calorimeter.FPCoolPipeTorRadius") ); //Tor radius of front plate outer cooling pipe
       calo_->G4Info_.set("FPCoolPipeRadius",       config.getDouble("calorimeter.FPCoolPipeRadius") );    //Inner radius of front plate outer cooling pipe
       calo_->G4Info_.set("FPCoolPipeThickness",    config.getDouble("calorimeter.FPCoolPipeThickness") ); //Thickness of front plate outer cooling pipe
       calo_->G4Info_.set("nPipes",                 config.getInt   ("calorimeter.nPipes") );              //Number of cooling pipes inside front plate
       calo_->G4Info_.set("pipeRadius",             config.getDouble("calorimeter.pipeRadius") );          //Inner radius of cooling pipes
       calo_->G4Info_.set("pipeThickness",          config.getDouble("calorimeter.pipeThickness") );       //Thickness of cooling pipes
       calo_->G4Info_.set("pipeInitSeparation",     config.getDouble("calorimeter.pipeInitSeparation") );  //Cooling pipe characteristics (see constructDiskCalorimeter)
       calo_->G4Info_.set("radSmTor",               config.getDouble("calorimeter.radSmTor") );            //Cooling pipe characteristics (see constructDiskCalorimeter)
       calo_->G4Info_.set("xsmall",                 config.getDouble("calorimeter.xsmall") );              //Cooling pipe characteristics (see constructDiskCalorimeter)
       calo_->G4Info_.set("xdistance",              config.getDouble("calorimeter.xdistance") );           //Cooling pipe characteristics (see constructDiskCalorimeter)
       calo_->G4Info_.set("rInnerManifold",         config.getDouble("calorimeter.rInnerManifold") );      //Cooling pipe characteristics (see constructDiskCalorimeter)
       calo_->G4Info_.set("pipeTorRadius",          getVDouble(config,"pipeTorRadius", calo_->G4Info_.get<int>("nPipes")) ); //Tor radius of cooling pipes and other
       calo_->G4Info_.set("largeTorPhi",            getVDouble(config,"largeTorPhi",   calo_->G4Info_.get<int>("nPipes")) ); //cooling pipe characteristics
       calo_->G4Info_.set("smallTorPhi",            getVDouble(config,"smallTorPhi",   calo_->G4Info_.get<int>("nPipes")) ); //(see constructDiskCalorimeter)
       calo_->G4Info_.set("yposition",              getVDouble(config,"yposition",     calo_->G4Info_.get<int>("nPipes")) );
       calo_->G4Info_.set("straightEndPhi",         getVDouble(config,"straightEndPhi",calo_->G4Info_.get<int>("nPipes")) );

       calo_->G4Info_.set("nCrates",                config.getInt   ("calorimeter.nCrates") );              //number of FEB crates per disk
       calo_->G4Info_.set("nBoards",                config.getInt   ("calorimeter.nBoards") );              //number of boards per crate
       calo_->G4Info_.set("FEBToDiskZOffset",       config.getDouble("calorimeter.FEBToDiskZOffset") );     //Z offset between front face of FEB and front face of disk
       calo_->G4Info_.set("crateXLength",           config.getDouble("calorimeter.crateXLength") );         //X length of crate
       calo_->G4Info_.set("crateYLength",           config.getDouble("calorimeter.crateYLength") );         //Y length of crate
       calo_->G4Info_.set("crateZLength",           config.getDouble("calorimeter.crateZLength") );         //Z length of crate
       calo_->G4Info_.set("crateFShieldThickness",  config.getDouble("calorimeter.crateFShieldThickness") );//Front shield thickness
       calo_->G4Info_.set("crateBShieldThickness",  config.getDouble("calorimeter.crateBShieldThickness") );//Bottom shield thickness
       calo_->G4Info_.set("crateBShieldLength",     config.getDouble("calorimeter.crateBShieldLength") );   //Bottom shield length
       calo_->G4Info_.set("crateTThickness",        config.getDouble("calorimeter.crateTThickness") );      //Thickness of top crate plate
       calo_->G4Info_.set("crateSThickness",        config.getDouble("calorimeter.crateSThickness") );      //Thickness of side crate plate
       calo_->G4Info_.set("crateFShieldYLength",    config.getDouble("calorimeter.crateFShieldYLength") );  //Y length front shield
       calo_->G4Info_.set("crateFShieldDeltaZ",     config.getDouble("calorimeter.crateFShieldDeltaZ") );   //Z difference betrween front shield and crate
       calo_->G4Info_.set("radiatorThickness",      config.getDouble("calorimeter.radiatorThickness") );    //radiator thickness
       calo_->G4Info_.set("radiatorZLength",        config.getDouble("calorimeter.radiatorZLength") );      //radiator Z length
       calo_->G4Info_.set("activeStripThickness",   config.getDouble("calorimeter.activeStripThickness") ); //active strip thickness
       calo_->G4Info_.set("passiveStripThickness",  config.getDouble("calorimeter.passiveStripThickness") );//passive strip thickness
       calo_->G4Info_.set("cratePhiAngles",         getVDouble(config,"cratePhiAngles",calo_->G4Info_.get<int>("nCrates")) ); //phi angle of the crates


       // CALORIMETER AND TRACKER ORIGIN
       double zTrackerCenter  =  config.getDouble("mu2e.detectorSystemZ0");
       double xOrigin         = -config.getDouble("mu2e.solenoidOffset");
       double zCaloStart      =  config.getDouble("calorimeter.caloMotherZ0");
       double zCaloEnd        =  config.getDouble("calorimeter.caloMotherZ1");
       calo_->G4Info_.set("caloOrigin", CLHEP::Hep3Vector(xOrigin,0,0.5*(zCaloStart+zCaloEnd)) ); //origin of the Calorimeter mother volume in Mu2e
       calo_->trackerCenter_ = CLHEP::Hep3Vector(xOrigin,0,zTrackerCenter);                   //origin of the Tracker mother volume in Mu2e


       makeIt();
       checkIt();
    }

    void DiskCalorimeterMaker::makeIt(void)
    {
        double vdThickness           = calo_->G4Info_.get<double>("vdThickness");
        double caloMotherZ0          = calo_->G4Info_.get<double>("caloMotherZ0");
        double caloMotherZ1          = calo_->G4Info_.get<double>("caloMotherZ1");
        double crystalHalfXY         = calo_->G4Info_.get<double>("crystalXYLength")/2.0;
        double wrapperHalfThick      = calo_->G4Info_.get<double>("wrapperThickness")/2.0;
        double crystalHalfZLength    = calo_->G4Info_.get<double>("crystalZLength")/2.0;
        double crystalCapHalfZLength = calo_->G4Info_.get<double>("crystalCapZLength")/2.0;
        double crystalCellRadius     = crystalHalfXY + 2.0*wrapperHalfThick;
        double innerDiskRadius       = calo_->G4Info_.get<double>("diskInAlRingRIn");
        double innerCrysRadius       = calo_->G4Info_.get<double>("diskCrystalRIn");
        double outerCrysRadius       = calo_->G4Info_.get<double>("diskCrystalROut");
        double outerCaseRadius       = calo_->G4Info_.get<double>("diskCaseRingROut");
        double diskOutRailROut       = calo_->G4Info_.get<double>("diskOutRailROut");
        auto   shimStepsInRowId      = calo_->G4Info_.get<std::vector<int>>("shimStepsInRowId");
        auto   shimStepsOutRowId     = calo_->G4Info_.get<std::vector<int>>("shimStepsOutRowId");
        double FPCarbonThick         = calo_->G4Info_.get<double>("FPCarbonZLength");
        double FPFoamThick           = calo_->G4Info_.get<double>("FPFoamZLength");
        double FPCoolPipeRadius      = calo_->G4Info_.get<double>("FPCoolPipeRadius");
        double FPpipeRadius          = calo_->G4Info_.get<double>("pipeRadius");
        double diskCaseHalfZLength   = crystalHalfZLength + crystalCapHalfZLength;
        double BPHoleHalfZ           = calo_->G4Info_.get<double>("BPHoleZLength")/2.0;
        double FEEBoxHalfZ           = calo_->G4Info_.get<double>("FEEZLength")/2.0;
        double FEEBoxThick           = calo_->G4Info_.get<double>("FEEBoxThickness");
        double BPPipeRadiusHigh      = calo_->G4Info_.get<double>("BPPipeRadiusHigh");
        double BPPipeHalfZOffset     = calo_->G4Info_.get<double>("BPPipeZOffset")/2.0;
        double crateZLength          = calo_->G4Info_.get<double>("crateZLength");
        double crateFShieldThick     = calo_->G4Info_.get<double>("crateFShieldThickness");
        double crateFShieldDeltaZ    = calo_->G4Info_.get<double>("crateFShieldDeltaZ");
        double FEBToDiskZOffset      = calo_->G4Info_.get<double>("FEBToDiskZOffset");
        double motherHalfZ           = (caloMotherZ1 - caloMotherZ0)/2.0;
        const auto& diskCrystalFile  = calo_->G4Info_.get<std::string>("diskCrystalFile");
        const auto& caloOrigin       = calo_->G4Info_.get<CLHEP::Hep3Vector>("caloOrigin");



        // First, calculate the total z length of the disk and the feb since we are not allowed
        // to get this from constructDiskCalorimeter.cc (no dependency on MC simulation)
        // Make sure this matches the geometry implemented in constructDiskCalorimeter !!!
        //
        double FPHalfZLength     = (FPCarbonThick + FPFoamThick - FPpipeRadius + FPCoolPipeRadius)/2.0;
        double BPHalfZLength     = BPHoleHalfZ + FEEBoxHalfZ + 2.0*FEEBoxThick + BPPipeHalfZOffset + BPPipeRadiusHigh;
        double diskHalfZLength   = FPHalfZLength + diskCaseHalfZLength + BPHalfZLength + vdThickness;
        double FEBZLength        = (crateZLength + crateFShieldDeltaZ + crateFShieldThick + 2.0*vdThickness);

        // Offset between front face of FEB and front face of disk
        double crateToDiskDeltaZ = FEBToDiskZOffset + vdThickness;

        // Offsets between the disk and crystal cordinate systems, i.e. distance between center of disk and front face crystals
        double disp = -diskHalfZLength + vdThickness + 2*FPHalfZLength + diskCaseHalfZLength - crystalHalfZLength + crystalCapHalfZLength;
        auto   diskOriginToCrystalOrigin = CLHEP::Hep3Vector(0,0,disp);


        // Then create the disks
        size_t crystalOffset(0);
        for (int idisk=0; idisk<CaloConst::_nDisk; ++idisk) {
           double separation    = calo_->G4Info_.get<std::vector<double>>("diskZMotherShift").at(idisk);
           double angleZ        = 0;

           double dR1 = innerDiskRadius;
           double dR2 = diskOutRailROut;
           double dZ  = 2.0*diskHalfZLength;

           CLHEP::Hep3Vector size(dR1,dR2,dZ) ;
           CLHEP::Hep3Vector originLocal(0, 0, -motherHalfZ + diskHalfZLength + separation + crateToDiskDeltaZ);

           CLHEP::Hep3Vector frontFaceCenter = caloOrigin + originLocal + diskOriginToCrystalOrigin;
           CLHEP::Hep3Vector backFaceCenter  = frontFaceCenter + CLHEP::Hep3Vector(0,0,2.0*crystalHalfZLength);
           CLHEP::HepRotation diskRotation   = CLHEP::HepRotation::IDENTITY*CLHEP::HepRotationZ(angleZ);

           auto thisDisk = Disk(idisk,innerCrysRadius,outerCrysRadius, 2.0*crystalCellRadius,
                                2.0*crystalHalfZLength, crystalOffset, diskOriginToCrystalOrigin,
                                diskCrystalFile);

           thisDisk.diskInfo().originLocal(originLocal);
           thisDisk.diskInfo().setPose(caloOrigin + originLocal, diskRotation);
           thisDisk.diskInfo().frontFaceCenter(frontFaceCenter);
           thisDisk.diskInfo().backFaceCenter(backFaceCenter);
           thisDisk.diskInfo().crystalDirection(CLHEP::Hep3Vector(0,0,1));
           thisDisk.diskInfo().size(size);
           thisDisk.diskInfo().originToCrystalOrigin(diskOriginToCrystalOrigin);
           thisDisk.diskInfo().FEBZOffset(crateToDiskDeltaZ);
           thisDisk.diskInfo().FEBZLength(FEBZLength);
           thisDisk.diskInfo().envelopeRad(dR1,dR2);
           thisDisk.diskInfo().crystalZlength(2.0*crystalHalfZLength);

           //fill the full Crystal List / diskId (direct access for performance optimization)
           for (unsigned icry=0;icry<thisDisk.nCrystals();++icry){
              Crystal& thisCrystal = thisDisk.crystal(icry);
              calo_->crystals_.push_back(&thisCrystal);

              //precompute the neighbors in the global frame
              thisCrystal.setID(icry+crystalOffset);
              thisCrystal.setNeighbors(thisDisk.neighbors(icry,1));
              thisCrystal.setNextNeighbors(thisDisk.neighbors(icry,2));

              //pre-compute the crystal position in the mu2e (global) frame
              CLHEP::Hep3Vector globalPosition = thisDisk.diskInfo().origin() + thisDisk.diskInfo().inverseRotation()*(thisCrystal.localPosition());
              thisCrystal.setPosition(globalPosition);
           }

           //calculate the position of the inner and outer steps, including the special shims to make the walls straight
           std::vector<double> par,stepsInX, stepsInY, stepsOutX, stepsOutY;
           int nrows = int(outerCaseRadius/2.0/crystalCellRadius)+1;
           for (int i=-nrows;i<nrows;++i) {

             auto cryList = thisDisk.idxFromRow(i);
             if (cryList.empty()) continue;

             auto sortFunctor = [&thisDisk](const auto i, const auto j)
                                   {return thisDisk.crystal(i).localPosition().x() <
                                           thisDisk.crystal(j).localPosition().x();};
             std::sort(std::begin(cryList),std::end(cryList),sortFunctor);

             // calculate the ymin,ymax, xmin/xmax for positive and negative side
             // there is a gap if the difference between the max negative and min positive > crystal size
             double ymin(1e6),ymax(-1e6),xLowMin(1e6),xLowMax(-1e6),xUpMin(1e6),xUpMax(-1e6);
             for (size_t ic=0;ic<cryList.size();++ic) {
               const auto& crystalI = thisDisk.crystal(cryList[ic]);
               ymin = std::min(ymin,crystalI.localPosition().y()-crystalI.size().y()/2.0);
               ymax = std::max(ymax,crystalI.localPosition().y()+crystalI.size().y()/2.0);
               if (crystalI.localPosition().x()<0.0){
                 xLowMin = std::min(xLowMin,crystalI.localPosition().x()-crystalI.size().x()/2.0);
                 xLowMax = std::max(xLowMax,crystalI.localPosition().x()+crystalI.size().x()/2.0);
               } else{
                 xUpMin = std::min(xUpMin,crystalI.localPosition().x()-crystalI.size().x()/2.0);
                 xUpMax = std::max(xUpMax,crystalI.localPosition().x()+crystalI.size().x()/2.0);
               }
             }

             //Now there are extra shims in the model, must add them
             double p0(xLowMin),p1(xUpMax);
             if (std::find(shimStepsOutRowId.begin(),shimStepsOutRowId.end(),i) != shimStepsOutRowId.end()){
               p0 -= crystalCellRadius;
               p1 += crystalCellRadius;
             }
             stepsOutX.push_back(p0);
             stepsOutX.push_back(p1);
             stepsOutY.push_back(ymin);
             stepsOutY.push_back(ymax);

             if (xUpMin - xLowMax > crystalCellRadius) {
                double p4(xLowMax),p5(xUpMin);
                if (std::find(shimStepsInRowId.begin(),shimStepsInRowId.end(),i) != shimStepsInRowId.end()){
                  p4 += crystalCellRadius;
                  p5 -= crystalCellRadius;
                }
                stepsInX.push_back(p4);stepsInX.push_back(p5);
                stepsInY.push_back(ymin);stepsInY.push_back(ymax);
             }
           }
           calo_->G4Info_.set("stepsInsideX",stepsInX);
           calo_->G4Info_.set("stepsInsideY",stepsInY);
           calo_->G4Info_.set("stepsOutsideX",stepsOutX);
           calo_->G4Info_.set("stepsOutsideY",stepsOutY);

           if (verbosityLevel_) std::cout<<"Constructed Disk "<<thisDisk.id()<<":  Rin="<<thisDisk.diskInfo().innerEnvelopeR()
                                         <<"  Rout="<<thisDisk.diskInfo().outerEnvelopeR()
                                         <<" (X,Y,Z)="<<thisDisk.diskInfo().origin()<<"  local_(X,Y,Z)="<<thisDisk.diskInfo().originLocal()
                                         <<"  with "<<thisDisk.nCrystals()<<" crystals"<<std::endl;

           if (verbosityLevel_ > 1) thisDisk.print();

           crystalOffset += thisDisk.nCrystals();
           calo_->disks_.push_back(std::move(thisDisk));
        }
    }



    void DiskCalorimeterMaker::checkIt(void)
    {

        //Consistency check between disk and CaloConst content
        for (size_t idisk=0;idisk<calo_->nDisks();++idisk)
        {
           const auto& thisDisk = calo_->disk(idisk);
           if (thisDisk.nCrystals() != CaloConst::_nCrystalPerDisk)
           throw cet::exception("DiskCaloGeom") << "The number of crystals ("<<thisDisk.nCrystals()
                                                <<")is differfent from "<<CaloConst::_nCrystalPerDisk<<"\n";
        }

        //check calorimeter fits inside mother envelope
        double diskRin       = calo_->G4Info_.get<double>("diskInAlRingRIn");
        double diskRout      = calo_->G4Info_.get<double>("diskCaseRingROut");
        double outerRingROut = calo_->G4Info_.get<double>("diskOutRailROut");
        if (diskRin > diskRout)
              throw cet::exception("DiskCaloGeom") << "calorimeter.diskInnerRingIn > calorimeter.diskOuterRingOut \n";

        if (outerRingROut > calo_->G4Info_.get<double>("caloDiskRadiusOut"))
              throw cet::exception("DiskCaloGeom") << "calorimeter outer radius larger than calorimeter mother \n";

        if (diskRin < calo_->G4Info_.get<double>("caloDiskRadiusIn"))
              throw cet::exception("DiskCaloGeom") << "calorimeter inner radius smaller than calorimeter mother \n";


        //check that holes in back plate are smaller than crystal, RO smaller than holes and FEE boxes fit
        if (calo_->G4Info_.get<double>("BPHoleXLength") > calo_->G4Info_.get<double>("crystalXYLength") ||
            calo_->G4Info_.get<double>("BPHoleYLength") > calo_->G4Info_.get<double>("crystalXYLength") )
              throw cet::exception("DiskCaloGeom") << "calorimeter backplate hole greater than crystal dimensions in X or Y \n";

        if (calo_->G4Info_.get<double>("readoutXLength") > calo_->G4Info_.get<double>("BPHoleXLength")  ||
            calo_->G4Info_.get<double>("readoutYLength") > calo_->G4Info_.get<double>("BPHoleYLength"))
              throw cet::exception("DiskCaloGeom") << "calorimeter readout larger than hole in X or Y \n";

        if (calo_->G4Info_.get<double>("readoutZLength") > calo_->G4Info_.get<double>("BPHoleZLength"))
              throw cet::exception("DiskCaloGeom") << "calorimeter readout too thick to fit in hole \n";

        if (calo_->G4Info_.get<double>("FEEXLength") > calo_->G4Info_.get<double>("BPHoleXLength")/calo_->G4Info_.get<int>("nSiPMPerCrystal"))
              throw cet::exception("DiskCaloGeom") << "calorimeter FEE box does not fit in X direction \n";

        if (calo_->G4Info_.get<double>("FEEYLength") > calo_->G4Info_.get<double>("crystalXYLength")-2*calo_->G4Info_.get<double>("FEEBoxThickness"))
              throw cet::exception("DiskCaloGeom") << "calorimeter FEE box does not fit in Y direction \n";


        //Check pipes
        for (int i=0;i<calo_->G4Info().get<int>("nPipes");++i){
           double pipeTorRadius     = calo_->G4Info_.get<std::vector<double>>("pipeTorRadius").at(i);
           double pipeRadius        = calo_->G4Info_.get<double>("pipeRadius");
           double radiusIn          = calo_->G4Info_.get<double>("FPInnerRadius");
           double radiusOut         = calo_->G4Info_.get<double>("FPOuterRadius");
           double foamZLength       = calo_->G4Info_.get<double>("FPFoamZLength");
           double carbonThick       = calo_->G4Info_.get<double>("FPCarbonZLength");
           double foamThick         = calo_->G4Info_.get<double>("FPFoamZLength");
           double coolPipeTorRadius = calo_->G4Info_.get<double>("FPCoolPipeTorRadius");
           double coolPipeRadius    = calo_->G4Info_.get<double>("FPCoolPipeRadius");

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
        int nBoards           = calo_->G4Info_.get<int>("nBoards");
        double radiatorDY     = calo_->G4Info_.get<double>("radiatorThickness")/2.0;
        double activeStripDY  = calo_->G4Info_.get<double>("activeStripThickness")/2.0;
        double passiveStripDY = calo_->G4Info_.get<double>("passiveStripThickness")/2.0;
        double crateXLength   = calo_->G4Info_.get<double>("crateXLength");
        double crateYLength   = calo_->G4Info_.get<double>("crateYLength");

        if ( nBoards*(radiatorDY+activeStripDY+passiveStripDY) > crateYLength)
              throw cet::exception("DiskCaloGeom") << "calorimeter FEB boards too thick\n";

        if (calo_->G4Info_.get<double>("crateFShieldYLength") > crateYLength)
              throw cet::exception("DiskCaloGeom") << "calorimeter FEB front shile too long in Y direction\n";

        if (calo_->G4Info_.get<double>("crateSThickness") > crateXLength)
              throw cet::exception("DiskCaloGeom") << "calorimeter FEB crate side too thick\n";

        if (calo_->G4Info_.get<double>("crateTThickness") > crateYLength )
              throw cet::exception("DiskCaloGeom") << "calorimeter FEB crate top too thick\n";

    }


    std::vector<double> DiskCalorimeterMaker::getVDouble(const SimpleConfig& config, const std::string& key, int size)
    {
       std::vector<double> vec;
       if (size>0 ) config.getVectorDouble("calorimeter."+key, vec, size);
       else         config.getVectorDouble("calorimeter."+key, vec);
       return vec;
    }

    std::vector<int> DiskCalorimeterMaker::getVInt(const SimpleConfig& config, const std::string& key, int size)
    {
       std::vector<int> vec;
       if (size>0) config.getVectorInt("calorimeter."+key, vec, size);
       else        config.getVectorInt("calorimeter."+key, vec);
       return vec;
    }


}
