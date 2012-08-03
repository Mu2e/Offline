// Andrei Gaponenko, 2011

#include "ExtinctionMonitorFNAL/inc/ExtMonFNAL_Maker.hh"
#include "ExtinctionMonitorFNAL/inc/ExtMonFNALMagnetMaker.hh"

#include "ExtinctionMonitorFNAL/inc/ExtMonFNAL.hh"
#include "ExtinctionMonitorFNAL/inc/ExtMonFNALBuilding.hh"

#include "cetlib/exception.h"
#include "boost/range/algorithm_ext/is_sorted.hpp"

#include "CLHEP/Units/SystemOfUnits.h"

#include "ConfigTools/inc/SimpleConfig.hh"

//#define AGDEBUG(stuff) std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<": "<<stuff<<std::endl;
#define AGDEBUG(stuff)

namespace mu2e {
  namespace ExtMonFNAL {

    //================================================================
    ExtMonFNALSensorStack ExtMonMaker::readStack(const SimpleConfig& config,
                                                 const std::string& prefix,
                                                 const CLHEP::Hep3Vector& refPointInMu2e,
                                                 const CLHEP::HepRotation& rotationInMu2e
                                                 )
    {
      ExtMonFNALSensorStack st;
      st.m_stackRefPointInMu2e = refPointInMu2e;
      st.m_stackRotationInMu2e = rotationInMu2e;
      st.m_coordinateRotationInMu2e = st.m_stackRotationInMu2e.inverse();

      config.getVectorDouble(prefix+".sensor_zoffset", st.m_sensor_zoffset, -1);
      config.getVectorDouble(prefix+".sensor_xoffset", st.m_sensor_xoffset, st.m_sensor_zoffset.size());
      config.getVectorDouble(prefix+".sensor_yoffset", st.m_sensor_yoffset, st.m_sensor_zoffset.size());
      config.getVectorDouble(prefix+".sensor_halfdx",  st.m_sensor_halfdx,  st.m_sensor_zoffset.size());
      config.getVectorDouble(prefix+".sensor_halfdy",  st.m_sensor_halfdy,  st.m_sensor_zoffset.size());
      config.getVectorDouble(prefix+".sensor_halfdz",  st.m_sensor_halfdz,  st.m_sensor_zoffset.size());
      config.getVectorDouble(prefix+".readout_halfdz", st.m_readout_halfdz, st.m_sensor_zoffset.size());

      if(!boost::is_sorted(st.m_sensor_zoffset)) {
        throw cet::exception("GEOM")<<"ExtMonFNAL_Maker: ERROR: "
                                    <<prefix<<".sensor_zoffset must be sorted in the ascending order"
                                    <<"\n"
          ;
      }

      //----------------------------------------------------------------
      // Test material plates

      config.getVectorString(prefix+".testMaterial.materials", st.m_testMaterialNames);
      if(!st.m_testMaterialNames.empty()) {
        const double testMaterialThickness = config.getDouble(prefix+".testMaterial.thickness");
        const double testMaterialXY = config.getDouble(prefix+".testMaterial.halfsize");

        st.m_testMaterialHalfSize.resize(3);
        st.m_testMaterialHalfSize[0] = testMaterialXY;
        st.m_testMaterialHalfSize[1] = testMaterialXY;
        st.m_testMaterialHalfSize[2] = testMaterialThickness/2;

        st.distanceToTestMaterials_ = config.getDouble(prefix+".testMaterial.distanceToStack");
        st.m_testMaterialPitch = testMaterialThickness + config.getDouble(prefix+".testMaterial.spacing");
      }

      //----------------------------------------------------------------
      return st;
    }

    //================================================================
    std::auto_ptr<ExtMon> ExtMonMaker::make(const SimpleConfig& config, const ExtMonFNALBuilding& room)
    {
      AGDEBUG("ExtMonFNAL make() begin");
      const int verbose = config.getInt("extMonFNAL.verbosityLevel", 0);

      std::auto_ptr<ExtMon> det(new ExtMon());

      //----------------------------------------------------------------
      // The upstream stack

      const double detectorDistanceToWall = config.getDouble("extMonFNAL.detectorDistanceToWall");
      const CLHEP::Hep3Vector upRefPointInMu2e = room.filterExitInMu2e()
        + room.collimator2RotationInMu2e() * CLHEP::Hep3Vector(0, 0, -detectorDistanceToWall);

      det->up_ = readStack(config, "extMonFNAL.up", upRefPointInMu2e, room.collimator2RotationInMu2e());

      //----------------------------------------------------------------
      // Spectrometer magnet

      const double upToMagnet = config.getDouble("extMonFNAL.up.distanceToMagnet");
      const CLHEP::Hep3Vector magnetRefInMu2e = det->up_.m_stackRefPointInMu2e
        + det->up_.m_stackRotationInMu2e * CLHEP::Hep3Vector(0,0, -upToMagnet);

      const double dp = config.getDouble("extMonFNAL.spectrometer.nominalMomentumAdjustment");
      det->spectrometerMagnet_ = ExtMonFNALMagnetMaker::read(config,
                                                             "extMonFNAL.spectrometer.magnet",
                                                             magnetRefInMu2e,
                                                             room.collimator2RotationInMu2e(),
                                                             room.filterMagnet().nominalMomentum() + dp
                                                             );

      //----------------------------------------------------------------
      // The downstream stack

      const double dnToMagnet = config.getDouble("extMonFNAL.dn.distanceToMagnet");
      const CLHEP::Hep3Vector dnRefPointInMu2e = det->spectrometerMagnet_.refPointInMu2e()
        + det->spectrometerMagnet().outRotationInMu2e() * CLHEP::Hep3Vector(0, 0, -dnToMagnet);

      det->dn_ = readStack(config, "extMonFNAL.dn", dnRefPointInMu2e, det->spectrometerMagnet_.outRotationInMu2e());

      //----------------------------------------------------------------
      if(verbose) {
        std::cout<<"ExtMonFNAL_Maker: UP stack center in Mu2e = "<<det->up_.m_stackRefPointInMu2e<<std::endl;
        std::cout<<"ExtMonFNAL_Maker: UP stackRotationInMu2e = "<<det->up_.m_stackRotationInMu2e<<std::endl;

        std::cout<<"ExtMonFNAL_Maker: magnet ref in Mu2e = "<<det->spectrometerMagnet_.refPointInMu2e()<<std::endl;
        std::cout<<"ExtMonFNAL_Maker: magnet center in Mu2e = "<<det->spectrometerMagnet_.geometricCenterInMu2e()<<std::endl;
        std::cout<<"ExtMonFNAL_Maker: magnet rotation in Mu2e = "<<det->spectrometerMagnet_.magnetRotationInMu2e()<<std::endl;

        std::cout<<"ExtMonFNAL_Maker: DN stack center in Mu2e = "<<det->dn_.m_stackRefPointInMu2e<<std::endl;
        std::cout<<"ExtMonFNAL_Maker: DN stackRotationInMu2e = "<<det->dn_.m_stackRotationInMu2e<<std::endl;
      }

      AGDEBUG("ExtMonFNAL maker end");

      return det;

    } // make()

  } // namespace ExtMonFNAL
} // namespace mu2e
