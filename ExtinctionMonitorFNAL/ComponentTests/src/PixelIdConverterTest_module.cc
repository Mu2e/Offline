// Andrei Gaponenko, 2012

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALPixelIdConverter.hh"
#include "DataProducts/inc/ExtMonFNALPixelId.hh"
#include "DataProducts/inc/ExtMonFNALPixelDenseId.hh"

#include <iostream>
#include <cstdlib>

#include "cetlib_except/exception.h"

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"

#include "CLHEP/Random/RandFlat.h"

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "SeedService/inc/SeedService.hh"

namespace mu2e {

  //================================================================
  class PixelIdConverterTest : public art::EDAnalyzer {
    unsigned int numTries_;
    CLHEP::RandFlat randFlat_;

    void test(const ExtMonFNALPixelIdConverter &conv, unsigned pix);

    void testException(const ExtMonFNALPixelIdConverter &conv, unsigned pix);
    void testException(const ExtMonFNALPixelIdConverter &conv, const ExtMonFNALPixelId& id);

  public:
    explicit PixelIdConverterTest(const fhicl::ParameterSet& pset);
    virtual void analyze(const art::Event& event);
  };

  //================================================================
  PixelIdConverterTest::PixelIdConverterTest(const fhicl::ParameterSet& pset)
    : art::EDAnalyzer(pset),
    , numTries_(pset.get<unsigned>("numTries"))
    , randFlat_(createEngine(art::ServiceHandle<SeedService>()->getSeed()))
  {}

  //================================================================
  void PixelIdConverterTest::analyze(const art::Event&) {

    GeomHandle<ExtMonFNAL::ExtMon> extmon;

    ExtMonFNALPixelIdConverter conv(*extmon);

    const unsigned totalNumberOfPixels = conv.totalNumberOfPixels();

    std::cout<<"PixelIdConverterTest: totalNumberOfPixels = "<<totalNumberOfPixels<<std::endl;

    std::cout<<"PixelIdConverterTest: testing edge cases"<<std::endl;
    test(conv, 0);
    test(conv, totalNumberOfPixels-1);

    std::cout<<"PixelIdConverterTest: performing = "<<numTries_<<" conversions"<<std::endl;
    for(unsigned i=0; i<numTries_; ++i) {
      test(conv, totalNumberOfPixels *  randFlat_.fire());
    }

    std::cout<<"ExtMonFNALPixelIdConverter: all tests PASSED"<<std::endl;
  }

  //================================================================
  void PixelIdConverterTest::test(const ExtMonFNALPixelIdConverter& conv, unsigned ipix) {

    ExtMonFNALPixelDenseId pix(ipix);
    ExtMonFNALPixelId id(conv.pixelId(pix));
    ExtMonFNALPixelDenseId pix2(conv.densePixelNumber(id));

    // std::cout<<pix<<", "<<id<<", "<<pix2<<std::endl;

    if(pix != pix2) {
      std::cerr<<"PixelIdConverterTest: ERROR: pix="<<pix<<", "<<id<<", pix2="<<pix2<<std::endl;
      abort();
    }
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::PixelIdConverterTest);
