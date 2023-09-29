//
//  Define the CRV test module geometry for KinKal
//  original author: David Brown (LBNL) 2023
//
#ifndef KinKalGeom_TestCRV_hh
#define KinKalGeom_TestCRV_hh
#include "KinKal/Geometry/Rectangle.hh"
#include <vector>
#include <memory>
namespace mu2e {
  namespace KinKalGeom {
    class TestCRV {
      public:
        using RecPtr = std::shared_ptr<KinKal::Rectangle>;

        // default constructor with nominal geometry
        TestCRV();
        // accessors
        // return by reference
        auto const& ex1() const { return *ex1_; }
        auto const& t1() const { return *t1_; }
        auto const& t2() const { return *t2_; }
        // return by ptr
        auto const& ex1Ptr() const { return ex1_; }
        auto const& t1Ptr() const { return t1_; }
        auto const& t2Ptr() const { return t2_; }

      private:
        RecPtr ex1_, t1_, t2_;
    };
  }
}

#endif
