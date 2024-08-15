// Andrei Gaponenko, 2011

#include "Offline/GeometryService/inc/ExtMonFNAL_Maker.hh"
#include "Offline/GeometryService/inc/ExtMonFNALMagnetMaker.hh"

#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALBuilding.hh"

#include "cetlib_except/exception.h"
#include "boost/range/algorithm_ext/is_sorted.hpp"

#include "CLHEP/Units/SystemOfUnits.h"

#include "Offline/ConfigTools/inc/SimpleConfig.hh"

//#define AGDEBUG(stuff) std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<": "<<stuff<<std::endl;
#define AGDEBUG(stuff)

namespace mu2e {
  namespace ExtMonFNAL {

    //================================================================
    ExtMonFNALPlaneStack ExtMonMaker::readStack(const SimpleConfig& config,
                                                const std::string& prefix,
                                                const CLHEP::Hep3Vector& refPointInMu2e,
                                                const CLHEP::HepRotation& rotationInMu2e,
                                                const ExtMonFNALModule& module
                                                )
    {
      ExtMonFNALPlaneStack pt;
      pt.m_stackRefPointInMu2e = refPointInMu2e;
      pt.m_stackRotationInMu2e = rotationInMu2e;
      pt.m_coordinateRotationInMu2e = pt.m_stackRotationInMu2e.inverse();
      config.getVectorDouble(prefix+".plane_zoffset", pt.m_plane_zoffset, -1);
      config.getVectorDouble(prefix+".plane_xoffset", pt.m_plane_xoffset, -1);
      config.getVectorDouble(prefix+".plane_yoffset", pt.m_plane_yoffset, -1);

      if(!boost::is_sorted(pt.m_plane_zoffset)) {
        throw cet::exception("GEOM")<<"ExtMonFNAL_Maker: ERROR: "
                                    <<prefix<<".plane_zoffset must be sorted in the ascending order"
                                    <<"\n";
      }

      std::vector<double> hs;
      config.getVectorDouble("extMonFNAL.planeHalfSize",  hs, 3);
      for(unsigned iplane = 0; iplane < pt.m_plane_zoffset.size(); ++iplane)
        {
          pt.planes_.emplace_back(module, hs);
          config.getVectorDouble(prefix+".module_zoffset", pt.planes_[iplane].m_module_zoffset, -1);
          config.getVectorDouble(prefix+".module_xoffset", pt.planes_[iplane].m_module_xoffset, -1);
          config.getVectorDouble(prefix+".module_yoffset", pt.planes_[iplane].m_module_yoffset, -1);
          config.getVectorDouble(prefix+".module_rotation", pt.planes_[iplane].m_module_rotation, -1);
          for(unsigned imodule = 0; imodule < pt.nModulesPerPlane(); imodule++)
            pt.planes_[iplane].m_module_rotation[imodule] *= CLHEP::degree;
        }

      return pt;
    }

    //================================================================
    std::unique_ptr<ExtMon> ExtMonMaker::make(const SimpleConfig& config, const ExtMonFNALBuilding& room)
    {
      AGDEBUG("ExtMonFNAL make() begin");
      const int verbose = config.getInt("extMonFNAL.verbosityLevel", 0);

      std::unique_ptr<ExtMon> det(new ExtMon());

      //----------------------------------------------------------------
      // Spectrometer magnet

      const double magnetToCollimatorDistance = config.getDouble("extMonFNAL.spectrometer.magnet.distanceToUpstreamWall")
        / (cos(room.filterAngleH())*cos(room.collimator2().angleV()));

      const CLHEP::Hep3Vector refEntranceToMagnet = room.filterExitInMu2e()
        + room.collimator2RotationInMu2e() * CLHEP::Hep3Vector(0, 0, -magnetToCollimatorDistance);

      const double dp = config.getDouble("extMonFNAL.spectrometer.nominalMomentumAdjustment");
      det->spectrometerMagnet_ = ExtMonFNALMagnetMaker::read(config,
                                                             "extMonFNAL.spectrometer.magnet",
                                                             room.collimator2RotationInMu2e(),
                                                             refEntranceToMagnet,
                                                             room.filterMagnet().nominalMomentum() + dp
                                                             );

      det->dnToExtMonCoordinateRotation_ =
        CLHEP::HepRotationX( -2 * det->spectrometerMagnet_.nominalBendHalfAngle());

      //----------------------------------------------------------------
      // Detector Mother
      config.getVectorDouble("extMonFNAL.detectorMotherHS", det->detectorMotherHS_, -1);

      double detectorMotherDistToMagnet = config.getDouble("extMonFNAL.detectorMotherDistToMagnet");
      double detectorMotherZCoord = det->detectorMotherHS()[1] - detectorMotherDistToMagnet - det->spectrometerMagnet().outerHalfSize()[1];
      CLHEP::Hep3Vector detectorMotherZVec = det->spectrometerMagnet().magnetRotationInMu2e()*CLHEP::Hep3Vector(0, detectorMotherZCoord, 0);
      det->detectorMotherCenterInMu2e_ = det->spectrometerMagnet().geometricCenterInMu2e() + detectorMotherZVec;


      //----------------

      config.getVectorDouble("extMonFNAL.sensorHalfSize", det->module_.sensorHalfSize_, -1);
      config.getVectorDouble("extMonFNAL.chipHalfSize", det->module_.chipHalfSize_, -1);
      //----------------------------------------------------------------
      // The upstream stack

      const CLHEP::Hep3Vector upRefPointInMu2e = det->spectrometerMagnet_.refPointInMu2e();
      det->up_ = readStack(config, "extMonFNAL.up", upRefPointInMu2e, room.collimator2RotationInMu2e(), det->module_);

      //----------------------------------------------------------------
      // The downstream stack

      const CLHEP::Hep3Vector dnRefPointInMu2e =  det->spectrometerMagnet_.refPointInMu2e();
      det->dn_ = readStack(config, "extMonFNAL.dn", dnRefPointInMu2e, det->spectrometerMagnet_.outRotationInMu2e(), det->module_);

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
