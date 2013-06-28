// Andrei Gaponenko, 2011

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL_Maker.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALMagnetMaker.hh"

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALBuilding.hh"

#include "cetlib/exception.h"
#include "boost/range/algorithm_ext/is_sorted.hpp"

#include "CLHEP/Units/SystemOfUnits.h"

#include "ConfigTools/inc/SimpleConfig.hh"

//#define AGDEBUG(stuff) std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<": "<<stuff<<std::endl;
#define AGDEBUG(stuff)

namespace mu2e {
  namespace ExtMonFNAL {

    //================================================================
    ExtMonFNALPlaneStack ExtMonMaker::readStack(const SimpleConfig& config,
                                                 const std::string& prefix,
                                                 const CLHEP::Hep3Vector& refPointInMu2e,
                                                 const CLHEP::HepRotation& rotationInMu2e
                                                 )
    {
      ExtMonFNALPlaneStack pt;
      pt.m_stackRefPointInMu2e = refPointInMu2e;
      pt.m_stackRotationInMu2e = rotationInMu2e;
      pt.m_coordinateRotationInMu2e = pt.m_stackRotationInMu2e.inverse();

      config.getVectorDouble(prefix+".plane_zoffset", pt.m_plane_zoffset, -1);
      config.getVectorDouble(prefix+".plane_xoffset", pt.m_plane_xoffset, -1);
      config.getVectorDouble(prefix+".plane_yoffset", pt.m_plane_yoffset, -1);
      config.getVectorDouble(prefix+".module_zoffset", pt.m_module_zoffset, -1);
      config.getVectorDouble(prefix+".module_xoffset", pt.m_module_xoffset, -1);
      config.getVectorDouble(prefix+".module_yoffset", pt.m_module_yoffset, -1);
      config.getVectorDouble(prefix+".motherTransverseHalfSize", pt.m_motherTransverseHalfSize, -1);
      pt.m_motherStartZ = config.getDouble(prefix+".motherStartZ"); 
      pt.m_motherEndZ = config.getDouble(prefix+".motherEndZ");

      if(!boost::is_sorted(pt.m_plane_zoffset)) {
        throw cet::exception("GEOM")<<"ExtMonFNAL_Maker: ERROR: "
                                    <<prefix<<".plane_zoffset must be sorted in the ascending order"
                                    <<"\n";
      }

      //----------------------------------------------------------------
      return pt;
    }

    //================================================================
    std::unique_ptr<ExtMon> ExtMonMaker::make(const SimpleConfig& config, const ExtMonFNALBuilding& room)
    {
      AGDEBUG("ExtMonFNAL make() begin");
      const int verbose = config.getInt("extMonFNAL.verbosityLevel", 0);

      std::unique_ptr<ExtMon> det(new ExtMon());

      config.getVectorDouble("extMonFNAL.planeHalfSize",  det->plane_.halfSize_, 3);
      config.getVectorDouble("extMonFNAL.sensorHalfSize", det->module_.sensorHalfSize_, -1);
      config.getVectorDouble("extMonFNAL.chipHalfSize", det->module_.chipHalfSize_, -1);
      //----------------------------------------------------------------
      // The upstream stack

      const double detectorDistanceToWall = config.getDouble("extMonFNAL.detectorDistanceToWall");
      const CLHEP::Hep3Vector upRefPointInMu2e = room.filterExitInMu2e()
        + room.collimator2RotationInMu2e() * CLHEP::Hep3Vector(0, 0, -detectorDistanceToWall);

      det->up_ = readStack(config, "extMonFNAL.up", upRefPointInMu2e, room.collimator2RotationInMu2e());

      //----------------------------------------------------------------
      // Spectrometer magnet

      const CLHEP::Hep3Vector magnetRefInMu2e = upRefPointInMu2e;
      const double dp = config.getDouble("extMonFNAL.spectrometer.nominalMomentumAdjustment");
      det->spectrometerMagnet_ = ExtMonFNALMagnetMaker::read(config,
                                                             "extMonFNAL.spectrometer.magnet",
                                                             magnetRefInMu2e,
                                                             room.collimator2RotationInMu2e(),
                                                             room.filterMagnet().nominalMomentum() + dp
                                                             );

      det->dnToExtMonCoordinateRotation_ =
        CLHEP::HepRotationX( -2 * det->spectrometerMagnet_.nominalBendHalfAngle());

      //----------------------------------------------------------------
      // The downstream stack

      const CLHEP::Hep3Vector dnRefPointInMu2e = magnetRefInMu2e;
      det->dn_ = readStack(config, "extMonFNAL.dn", dnRefPointInMu2e, det->spectrometerMagnet_.outRotationInMu2e());

      //----------------------------------------------------------------
      det->dn_.planeNumberOffset_ = 0;
      det->up_.planeNumberOffset_ = det->dn_.nplanes();

      //----------------------------------------------------------------
      if(verbose) {
        std::cout<<"ExtMonFNAL_Maker: UP stack center in Mu2e = "<<det->up_.m_stackRefPointInMu2e<<std::endl;
        std::cout<<"ExtMonFNAL_Maker: UP stackRotationInMu2e = "<<det->up_.m_stackRotationInMu2e<<std::endl;

        std::cout<<"ExtMonFNAL_Maker: magnet ref in Mu2e = "<<det->spectrometerMagnet_.refPointInMu2e()<<std::endl;
        std::cout<<"ExtMonFNAL_Maker: magnet center in Mu2e = "<<det->spectrometerMagnet_.geometricCenterInMu2e()<<std::endl;

        std::cout<<"ExtMonFNAL_Maker: magnet half bend angle  = "
                 <<det->spectrometerMagnet_.trackBendHalfAngle(room.filterMagnet().nominalMomentum() + dp)<<std::endl;

        std::cout<<"ExtMonFNAL_Maker: magnet rotation in Mu2e = "<<det->spectrometerMagnet_.magnetRotationInMu2e()<<std::endl;

        std::cout<<"ExtMonFNAL_Maker: DN stack center in Mu2e = "<<det->dn_.m_stackRefPointInMu2e<<std::endl;
        std::cout<<"ExtMonFNAL_Maker: DN stackRotationInMu2e = "<<det->dn_.m_stackRotationInMu2e<<std::endl;
      }

      AGDEBUG("ExtMonFNAL maker end");

      return det;

    } // make()

  } // namespace ExtMonFNAL
} // namespace mu2e
