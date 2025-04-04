//
//  Define the CRV test module geometry for KinKal
//  original author: David Brown (LBNL) 2023
//
#ifndef KinKalGeom_CRV_hh
#define KinKalGeom_CRV_hh
#include "KinKal/Geometry/Rectangle.hh"
#include <vector>
#include <memory>
namespace mu2e {
  namespace KinKalGeom {
    class CRV {
    public:
      using RecPtr = std::shared_ptr<KinKal::Rectangle>;
      // default constructor with nominal geometry
      CRV();
      // accessors
      // return by reference
      auto const& r1() const { return *r1_; }
      auto const& r2() const { return *r2_; }
      auto const& r3() const { return *r3_; }
      auto const& r4() const { return *r4_; }
      auto const& r5() const { return *r5_; }
      auto const& r6() const { return *r6_; }
      auto const& l1() const { return *l1_; }
      auto const& l2() const { return *l2_; }
      auto const& l3() const { return *l3_; }
      auto const& t1() const { return *t1_; }
      auto const& t2() const { return *t2_; }
      auto const& t3() const { return *t3_; }
      auto const& t4() const { return *t4_; }
      auto const& t5() const { return *t5_; }
      auto const& e1() const { return *e1_; }
      auto const& e2() const { return *e2_; }
      auto const& u1() const { return *u1_; }
      auto const& d1() const { return *d1_; }
      auto const& d2() const { return *d2_; }
      auto const& d3() const { return *d3_; }
      auto const& d4() const { return *d4_; }
      auto const& c1() const { return *c1_; }
      auto const& c2() const { return *c2_; }
      auto const& c3() const { return *c3_; }
      auto const& c4() const { return *c4_; }
      // return by ptr
      auto const& r1Ptr() const { return r1_; }
      auto const& r2Ptr() const { return r2_; }
      auto const& r3Ptr() const { return r3_; }
      auto const& r4Ptr() const { return r4_; }
      auto const& r5Ptr() const { return r5_; }
      auto const& r6Ptr() const { return r6_; }
      auto const& l1Ptr() const { return l1_; }
      auto const& l2Ptr() const { return l2_; }
      auto const& l3Ptr() const { return l3_; }
      auto const& t1Ptr() const { return t1_; }
      auto const& t2Ptr() const { return t2_; }
      auto const& t3Ptr() const { return t3_; }
      auto const& t4Ptr() const { return t4_; }
      auto const& t5Ptr() const { return t5_; }
      auto const& e1Ptr() const { return e1_; }
      auto const& e2Ptr() const { return e2_; }
      auto const& u1Ptr() const { return u1_; }
      auto const& d1Ptr() const { return d1_; }
      auto const& d2Ptr() const { return d2_; }
      auto const& d3Ptr() const { return d3_; }
      auto const& d4Ptr() const { return d4_; }
      auto const& c1Ptr() const { return c1_; }
      auto const& c2Ptr() const { return c2_; }
      auto const& c3Ptr() const { return c3_; }
      auto const& c4Ptr() const { return c4_; }
      private:
      RecPtr r1_, r2_, r3_, r4_, r5_, r6_, l1_, l2_, l3_, t1_, t2_, t3_, t4_, t5_, e1_, e2_, u1_, d1_, d2_, d3_, d4_, c1_, c2_, c3_, c4_;
    };
  }
}

#endif
