//
//  Define the CRV test module geometry for KinKal
//  original author: David Brown (LBN) 2023
//
#ifndef KinKalGeom_TestCRV_hh
#define KinKalGeom_TestCRV_hh
#include "KinKal/Geometry/Rectangle.hh"
#include <vector>
namespace mu2e {
  namespace KinKalGeom {
    class TestCRV {
      public:
        // default constructor with nominal geometry
        TestCRV();
        // accessors
        auto const& modulesEX1() const { return ex1_; }
        auto const& modulesT1() const { return t1_; }
        auto const& modulesT2() const { return t2_; }

      private:
        KinKal::Rectangle ex1_, t1_, t2_;
    };
  }
}

#endif
