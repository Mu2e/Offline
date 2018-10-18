// Andrei Gaponenko, 2012

#include "ExtinctionMonitorFNAL/Reconstruction/inc/TrackExtrapolator.hh"
#include <iostream>
#include <cstdlib>

#include "cetlib_except/exception.h"

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"

#include "RecoDataProducts/inc/ExtMonFNALTrkParam.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "GeometryService/inc/GeomHandle.hh"

namespace mu2e {
  namespace ExtMonFNAL {

    //================================================================
    class TrackExtrapolatorTest : public art::EDAnalyzer {
      void testToPlane(const TrackExtrapolator& ex, const ExtMonFNALTrkParam& start, unsigned toPlane);
      void testThroughMagnet(const TrackExtrapolator& ex, const ExtMonFNALTrkParam& start);

    public:
      explicit TrackExtrapolatorTest(const fhicl::ParameterSet& pset);
      virtual void analyze(const art::Event& event);
    };

    //================================================================
    TrackExtrapolatorTest::TrackExtrapolatorTest(const fhicl::ParameterSet& pset)
      : art::EDAnalyzer(pset) {}

    //================================================================
    void TrackExtrapolatorTest::analyze(const art::Event&) {

      GeomHandle<ExtMonFNAL::ExtMon> extmon;

      TrackExtrapolator ex(&*extmon);
      const ExtMonFNALMagnet& mag = extmon->spectrometerMagnet();

      const double nominalHalfBendAngle = mag.trackBendHalfAngle(mag.nominalMomentum());

      std::cout<<"INFO: nominal bend half angle = "<<nominalHalfBendAngle<<std::endl;
      std::cout<<"INFO: nominal bend radius = "<<mag.trackBendRadius(mag.nominalMomentum())<<std::endl;
      std::cout<<"INFO: magnet half length = "<<mag.outerHalfSize()[2]<<std::endl;
      std::cout<<std::endl;

      std::cout<<"TrackExtrapolatorTest: extrapolating reference trajectory"<<std::endl;
      ExtMonFNALTrkParam ref;
      ref.setz0(extmon->up().sensor_zoffset().back());
      ref.setrinv(1./mag.trackBendRadius(mag.nominalMomentum()));
      testToPlane(ex, ref, 0);

      // Now raise the track by 20 mm
      std::cout<<"TrackExtrapolatorTest: extrapolating at y=+20"<<std::endl;
      ref.setposy(20.);
      testToPlane(ex, ref, 0);

      // point track to a side
      std::cout<<"TrackExtrapolatorTest: slopex"<<std::endl;
      ref.setslopex(0.005);
      testToPlane(ex, ref, 0);

      // and slightly down
      std::cout<<"TrackExtrapolatorTest: slopey down"<<std::endl;
      ref.setslopey(+0.5 * nominalHalfBendAngle);
      testToPlane(ex, ref, 0);

      // and up
      std::cout<<"TrackExtrapolatorTest: slopey up"<<std::endl;
      ref.setslopey(- 1. * nominalHalfBendAngle);
      testToPlane(ex, ref, 0);

      // adjust momentum so that it ~ comes back
      std::cout<<"TrackExtrapolatorTest: lower momentum by half to compensate for the slope"<<std::endl;
      ref.setrinv(2*ref.rinv());
      testToPlane(ex, ref, 0);

      // adjust momentum so that it ~ comes back
      std::cout<<"Same track in a single stack"<<std::endl;
      testToPlane(ex, ref, extmon->dn().nplanes());

      // Make track curl inside the magnet
      std::cout<<"TrackExtrapolatorTest: track that barely makes it through the magnet"<<std::endl;
      ExtMonFNALTrkParam curl;
      curl.setz0(extmon->up().sensor_zoffset()[0]);
      curl.setslopey(nominalHalfBendAngle);
      curl.setrinv(0.999/(2*mag.outerHalfSize()[2]));
      testThroughMagnet(ex, curl);

      // Stronger curl inside the magnet
      std::cout<<"TrackExtrapolatorTest: track that should not make it through the magnet"<<std::endl;
      curl.setrinv(1.001/(2*mag.outerHalfSize()[2]));
      testThroughMagnet(ex, curl);

    }

    //================================================================
    void TrackExtrapolatorTest::testToPlane(const TrackExtrapolator& ex, const ExtMonFNALTrkParam &start, unsigned toPlane) {

      std::cout<<"TrackExtrapolatorTest: start = "<<start<<std::endl;
      ExtMonFNALTrkParam pp(start);
      if(ex.extrapolateToPlane(toPlane, &pp)) {
        std::cout<<"TrackExtrapolatorTest:   end = "<<pp<<std::endl;
      }
      else {
        std::cout<<"TrackExtrapolatorTest: DID NOT GO through magnet, now at = "<<pp<<std::endl;
      }
    }

    //================================================================
    void TrackExtrapolatorTest::testThroughMagnet(const TrackExtrapolator& ex, const ExtMonFNALTrkParam &start) {

      ExtMonFNALTrkParam pp(start);
      std::cout<<"TestMagnet::     start = "<<pp<<std::endl;
      ex.extrapolateToMagnet(&pp);
      std::cout<<"TestMagnet:: at magnet = "<<pp<<std::endl;
      if(ex.extrapolateThroughMagnet(&pp)) {
        std::cout<<"TestMagnet::       end = "<<pp<<std::endl;
      }
      else {
        std::cout<<"TestMagnet:: DID NOT GO through magnet, now at = "<<pp<<std::endl;
      }
    }

    //================================================================
  } // namespace ExtMonFNAL
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::ExtMonFNAL::TrackExtrapolatorTest);
