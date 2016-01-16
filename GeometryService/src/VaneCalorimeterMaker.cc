//
// Make a Vane Calorimeter.
//
// original authors Julie Managan and Robert Bernstein

//
// C++ includes
#include <math.h>
#include <memory>


// Mu2e includes
#include "CalorimeterGeom/inc/VaneCalorimeter.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/Vane.hh"
#include "GeometryService/inc/VaneCalorimeterMaker.hh"

// Framework include files
#include "cetlib/exception.h"

// other includes
#include "CLHEP/Vector/RotationX.h"
#include "CLHEP/Vector/RotationY.h"
#include "CLHEP/Vector/RotationZ.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Matrix/DiagMatrix.h"


// Vane geometry
//
//    crystals + readout are places in the wrapping material
//    wrapper/shell are placed in the crystal box
//    shield / neutron absorber are in front of the crystal box
//    absorber + crystal box are placed inside the casing
//
//    cross section in z: casing - shield - absorber - wrapper - crystal - wrapper - crystal -... - wrapper - casing
//    cross section in height : casing - wrapper - crystal - readout - wrapper - casing
//    cross section in radius : casing - wrapper - crystal - wrapper - crystal - ...  - wrapper - casing


//
//
//    Local reference frame of the vane:
//      X is along the longitudinal section of the vane, i.e. the large size
//      Y is radial towards large radius
//      Z is along the crystal longitudinal side


//  Coordinate systems

//    Intro: We model the crtsyals as polyhedra to keep the generality with the disk case. In G4
//           polyhedra have the origin of the coordinate system at their _base_, and are always
//           oriented along the z direction.
//           On the other hand, the vanes have origin of the coordinate system at their _center_
//
//           This choice is the easiest for creation of the geometry and reconstruction.
//           It is also convenient when analyzing track/calo matching to have a coordinate system
//           at the front face of the active volume, the crystals in our case. This one is available
//           in addiiton to the base choice.

//           The idea is to go from the calorimeter -> vane frame -> crystal frame to navigate the geometry

//           __Here comes the annoying part. As mentioned, the crystals are oriented along the
//           Z axis. For matching tracks to clusters, it is convenient to have all type of
//           calorimeters with the z axis oriented along the crystal longitudinal size. This
//           means one wants a front view of the vane (i.e. seeing the matrix of crystals),
//           not a size view (i.e. facing the short side).

//           On the other hand, the z axis in the Mu2e coordinate frame would naturally be
//           defined as the axis along the large side of the vane. This means that the
//           vanes will have to be rotated to be properly aligned in Mu2e. Not a terrible
//           problem, except some quantities are defined in the Mu2e frame, while others in
//           the vane frame. So one need to pay attention when translating between frames! In
//           our case, x <-> z is the way of translating it!

//           One way to solve this is to still use a vane frame oriented along the Mu2e frame, the z
//           axis along the large side of the vanes, and apply an extra rotation to recover the vane frame
//           described above. The drawback is that the crystals must be rotated to fit into the vanes
//           (since their longitudinal sides are always oriented along the z coordinate system), and one need
//           an extra rotation for navigating from the vane -> mu2e frame.
//           If one use the "front face view" as the vane frame, this is a bit more annoying to code here
//           but better in other places, so bite the bullet here (I have indicated where this is messy)

//           In summary, the vane coordinate system is oriented in a "front face view", while the
//           "Mu2e view", aka"side view" is sometimes needed. Pay attention to this!



//    Coordinate system
//
//           - crystal:
//                      origin:      base of the crystal
//                      orientation: crystal oriented along the z direction (direction fron face - back face)
//                                   the front (back) face is at z=0 (z=crystal length), see above about polyhedra.
//           - vane:
//                      origin:      center of the vane
//                      orientation: x axis along larger side, y axis along smaller side, z along crystal longitudinal size (see face of the vane)
//                      other:       extra coordinate system placed at the front of the crystal - used for track-calo matching
//
//           - calorimeter:
//                      origin:      center of front face, point closest to the traclker along z
//
//         Note: the "vane" is defined as whatever volume containing the crystal is rotated in G4 (see DiskCalorimeterMaker for details)
//
//         Note: The extra disk coordinate system is suffixed "FF" for FrontFace
//
//    Placement:
//
//         The position of the crystals in the vanes (see vane.cc) = the position of the crystal w.r.t vane origin, given by
//         vaneOriginToCrystalOrigin:
//              vaneOriginToCrystalOrigin(-absorberHalfLength,0,-_calo->_caloGeomInfo.crystalHalfLength()-roHalfThickness);






namespace mu2e{

    VaneCalorimeterMaker::VaneCalorimeterMaker( SimpleConfig const& config, double solenoidOffset)
    {

        _calo = std::unique_ptr<VaneCalorimeter>(new VaneCalorimeter());
        _calo->_caloType = Calorimeter::CaloType::vane;

        _calo->_nSections             = config.getInt(   "calorimeter.numberOfVanes");
        _calo->_nCrystalX             = config.getInt   ("calorimeter.nCrystalXSlices");
        _calo->_nCrystalY             = config.getInt   ("calorimeter.nCrystalYSlices");
        _calo->_shieldHalfThickness   = config.getDouble("calorimeter.shieldHalfThickness");
        _calo->_absorberHalfThickness = config.getDouble("calorimeter.neutronAbsorberHalfThickness");

        //Fill the Common Calo Data
        _calo->_caloGeomInfo.nROPerCrystal(      config.getInt("calorimeter.crystalReadoutChannelCount"));
        _calo->_caloGeomInfo.crystalHalfTrans(   config.getDouble("calorimeter.crystalHalfTrans") );
        _calo->_caloGeomInfo.crystalHalfLength(  config.getDouble("calorimeter.crystalHalfLong") );
        _calo->_caloGeomInfo.wrapperThickness(   config.getDouble("calorimeter.crystalWrapperThickness") );
        _calo->_caloGeomInfo.roHalfTrans(        config.getDouble("calorimeter.crystalReadoutHalfTrans") );
        _calo->_caloGeomInfo.roHalfThickness(    config.getDouble("calorimeter.crystalReadoutHalfThickness") );
        _calo->_caloGeomInfo.caseThickness(      config.getDouble("calorimeter.caseThickness") );
        _calo->_caloGeomInfo.enveloppeInRadius(  config.getDouble("calorimeter.caloMotherInRadius") );
        _calo->_caloGeomInfo.enveloppeOutRadius( config.getDouble("calorimeter.caloMotherOutRadius") );
        _calo->_caloGeomInfo.enveloppeZ0(        config.getDouble("calorimeter.caloMotherZ0") );
        _calo->_caloGeomInfo.enveloppeZ1(        config.getDouble("calorimeter.caloMotherZ1") );

        _calo->_caloGeomInfo.apdMeanNoise(       config.getDouble("calorimeter.meanNoiseAPD", 0.0) );
        _calo->_caloGeomInfo.apdSigmaNoise(      config.getDouble("calorimeter.sigmaNoiseAPD", 0.03) );
        _calo->_caloGeomInfo.lysoLightYield(     config.getDouble("calorimeter.lysoLightYield", 2000.0) );
        _calo->_caloGeomInfo.apdQuantumEff(      config.getDouble("calorimeter.quantumEffAPD", 0.68) );
        _calo->_caloGeomInfo.apdCollectEff(      config.getDouble("calorimeter.lightCollectEffAPD", 0.11));
        _calo->_caloGeomInfo.nonUniformity(      config.getDouble("calorimeter.crystalNonUniformity",0.0) );
        _calo->_caloGeomInfo.timeGap(            config.getDouble("calorimeter.timeGap",100.0) );
        _calo->_caloGeomInfo.electronEdep(       config.getDouble("calorimeter.electronDepositionAPD",1000.0) );
        _calo->_caloGeomInfo.electronEmin(       config.getDouble("calorimeter.electronMinEnergyAPD",0.1) );

        _calo->_caloGeomInfo.nPipes(             config.getInt("calorimeter.nPipes",0));
        _calo->_caloGeomInfo.pipeRadius(         config.getDouble("calorimeter.pipeRadius",5) );
        _calo->_caloGeomInfo.pipeThickness(      config.getDouble("calorimeter.pipeThickness",0.5) );

        std::vector<double> temp;
        config.getVectorDouble("calorimeter.pipeTorRadius", temp, _calo->_caloGeomInfo.nPipes());
        _calo->_caloGeomInfo.pipeTorRadius( temp );


        double crystalFullTrans    = _calo->_caloGeomInfo.crystalHalfTrans() + _calo->_caloGeomInfo.wrapperThickness();
        _calo->_rMin               = config.getDouble("calorimeter.rInscribed");
        _calo->_rMax               = _calo->_rMin + 2.0*crystalFullTrans*_calo->_nCrystalY;

        _calo->_isVaneTilted       = config.getBool("calorimeter.tiltVane", 0);

        //THE CALORIMETER ORIGIN IS TAKEN AS THE POINT CLOSEST TO THE TRACKER IN MU2E COORDINATES
        double xOrigin             = -config.getDouble("mu2e.solenoidOffset");
        double zCaloFront          = config.getDouble("calorimeter.calorimeterZFront");
        double zTrackerCenter      = config.getDouble("mu2e.detectorSystemZ0");
        _calo->_origin             = CLHEP::Hep3Vector(xOrigin,0,zCaloFront);
        _calo->_trackerCenter      = CLHEP::Hep3Vector(xOrigin,0,zTrackerCenter);


         //standard formula to get the volume of the rectangle
        _calo->_caloGeomInfo.crystalVolume(8*_calo->_caloGeomInfo.crystalHalfTrans()*_calo->_caloGeomInfo.crystalHalfTrans()*_calo->_caloGeomInfo.crystalHalfLength());

        _verbosityLevel = config.getInt("calorimeter.verbosityLevel",0);

        //make sure above information is consistent
        CheckIt();

        // Create vanes
        MakeVanes();

      }


      VaneCalorimeterMaker::~VaneCalorimeterMaker() {}




      void VaneCalorimeterMaker::MakeVanes()
      {

          double crystalHalfLength  = _calo->_caloGeomInfo.crystalHalfLength();
          double crystalHalfTrans   = _calo->_caloGeomInfo.crystalHalfTrans();
          double roHalfThickness    = _calo->_caloGeomInfo.roHalfThickness();
          double caseThickness      = _calo->_caloGeomInfo.caseThickness();
          double wrapperThickness   = _calo->_caloGeomInfo.wrapperThickness();

          double absorberHalfLength =  _calo->_shieldHalfThickness + _calo->_absorberHalfThickness;
          double crystalFullTrans   =  crystalHalfTrans + wrapperThickness;


          double dX     = crystalFullTrans * _calo->_nCrystalX + absorberHalfLength + caseThickness;
          double dY     = crystalFullTrans * _calo->_nCrystalY + caseThickness;
          double dZ     = crystalHalfLength + wrapperThickness + roHalfThickness + caseThickness;

          double radius = _calo->_rMin + dY - caseThickness;
          double dphi   = 2.0*CLHEP::pi/float(_calo->_nSections);

          CLHEP::Hep3Vector size(dX,dY,dZ);

          // This is where the offsets between the different coordinate systems are set, see Note for full explanation
          // Seriously, read the note at the top before changing this! Really!
          CLHEP::Hep3Vector vaneOriginToCrystalOrigin(absorberHalfLength,0,-_calo->_caloGeomInfo.crystalHalfLength()-roHalfThickness);


          for (unsigned int i=0; i<_calo->_nSections; ++i )
          {
              double phi = float(i)*dphi;
              CLHEP::Hep3Vector originLocal(radius*cos(phi),radius*sin(phi),dX);

              std::shared_ptr<Vane> thisVane( new Vane(i,_calo->_rMin,_calo->_nCrystalX,_calo->_nCrystalY, size,
                                                       crystalFullTrans, crystalHalfLength, vaneOriginToCrystalOrigin) );
              _calo->_sections.push_back(thisVane);


              thisVane->setOriginLocal(      originLocal );
              thisVane->setOrigin(           originLocal + _calo->origin() );
              thisVane->setRotation(         (CLHEP::HepRotation::IDENTITY)*CLHEP::HepRotationY(CLHEP::pi/2)*CLHEP::HepRotationZ(-phi+CLHEP::pi/2));
              if (_calo->_isVaneTilted) thisVane->setRotation(CLHEP::HepRotation::IDENTITY);


	      CLHEP::Hep3Vector frontFaceCenter = _calo->origin() + originLocal +  thisVane->rotation()*CLHEP::Hep3Vector(-dX,0,0);
	      CLHEP::Hep3Vector backFaceCenter  = _calo->origin() + originLocal +  thisVane->rotation()*CLHEP::Hep3Vector(dX,0,0);
              thisVane->setFrontFaceCenter(frontFaceCenter);
              thisVane->setBackFaceCenter(backFaceCenter);
	      thisVane->setEnveloppeRad(_calo->_rMin, _calo->_rMax);


              //fill the full Crystal List (direct access to crystal from calorimeter as requested from users)
              int crystalOffset = _calo->_fullCrystalList.size();
              for (int icry=0;icry<thisVane->nCrystals();++icry)
              {
                 Crystal& thisCrystal = thisVane->crystal(icry);
                 _calo->_fullCrystalList.push_back(&thisCrystal);

                 //precompute the neighbors in the global frame
                 thisCrystal.setNeighbors(_calo->neighborsByLevel(icry+crystalOffset,1));
                 thisCrystal.setNextNeighbors(_calo->neighborsByLevel(icry+crystalOffset,2));
                 thisCrystal.setPosition(_calo->crystal(icry).position());

                 //calculate the crystal position in the mu2e frame (aka global frame), taken from BaseCalorimeter.cc
                 CLHEP::Hep3Vector globalPosition = thisVane->origin() + thisVane->inverseRotation()*(thisCrystal.localPosition());
                 thisCrystal.setPosition(globalPosition);

                 //std::cout<<"Crystal neighbors for cry="<<icry<<"  :";for (auto const& cry : thisCrystal.neighbors()) std::cout<<cry<<" ";std::cout<<std::endl;
                 //std::cout<<icry<<"  "<<thisVane->idxFromPosition(thisCrystal.localPosition().x(),thisCrystal.localPosition().y())<<std::endl;
              }

              if (_verbosityLevel) std::cout<<"Constructed Vane "<<thisVane->id()<<": Rin="<<_calo->_rMin<<"  Rout="<<_calo->_rMin+dY
                                            <<" dX="<<dX<<" dY = "<<dY<<"  dZ="<<dZ<<" (X,Y,Z)="<<thisVane->origin()<<"  local_(X,Y,Z)="<<thisVane->originLocal()
                                            <<"  with "<<thisVane->nCrystals()<<" crystals"<<std::endl;

              if (_verbosityLevel > 1) thisVane->print();
          }

      }


      void VaneCalorimeterMaker::CheckIt(void)
      {


            //check that calorimeter fits in the mother volume
            double crystalFullTrans   =  _calo->_caloGeomInfo.crystalHalfTrans() + _calo->_caloGeomInfo.wrapperThickness();
            double caseThickness      = _calo->_caloGeomInfo.caseThickness();
            double absorberHalfLength =  _calo->_shieldHalfThickness + _calo->_absorberHalfThickness;
            double dX                 = crystalFullTrans * _calo->_nCrystalX + absorberHalfLength + caseThickness;
            double dY                 = crystalFullTrans * _calo->_nCrystalY + caseThickness;
            
	    double calozBegin         = _calo->_origin.z();
            double calozEnd           = _calo->_origin.z() + 2*dX;


            if ( (_calo->_rMin + 2*dY) > _calo->_caloGeomInfo.enveloppeOutRadius())
                     {throw cet::exception("VaneCaloGeom") << "calorimeter outer radius larger than calorimeter mother \n";}

            if (  _calo->_rMin < _calo->_caloGeomInfo.enveloppeInRadius())
                     {throw cet::exception("VaneCaloGeom") << "calorimeter inner radius smaller than calorimeter mother \n";}

            if (calozBegin <  _calo->_caloGeomInfo.enveloppeZ0() || calozBegin >  _calo->_caloGeomInfo.enveloppeZ1())
                     {throw cet::exception("VaneCaloGeom") << "calorimeter.calorimeterZFront   outside calorimeter mother.\n";}

            if (calozEnd   > _calo->_caloGeomInfo.enveloppeZ1())
                     {throw cet::exception("VaneCaloGeom") << "calorimeter z-coordinate extends outside calorimeter mother.\n";}


            // Check number of readouts
            int nRO    = _calo->caloGeomInfo().nROPerCrystal();
            if( ! (nRO==1 || nRO==2 || nRO==4) )
              {throw cet::exception("VaneCaloGeom") << "calorimeter.crystalReadoutChannelCount can only be 1,2 or 4.\n";}


            // Check size of readouts
            if( nRO==1 ) {
               if( _calo->_caloGeomInfo.roHalfTrans() > _calo->_caloGeomInfo.crystalHalfTrans() )
                 {throw cet::exception("VaneCaloGeom") << "calorimeter.crystalReadoutHalfTrans > calorimeter.crystalHalfTrans.\n";}

            } else {
               if( _calo->_caloGeomInfo.roHalfTrans()> 0.5*_calo->_caloGeomInfo.crystalHalfTrans() )
                 {throw cet::exception("VaneCaloGeom") << "calorimeter.crystalReadoutHalfTrans > 0.5*calorimeter.crystalHalfTrans.\n";}

            }

      }


}

