//
//  Define the nominal calorimeter boundary and reference surfaces, used to extrapolate and sample KinKal track fits, and to build
//  the passive materials in the fit
//  original author: Sophie Middleton (2025)
//
#ifndef KinKalGeom_Calo_hh
#define KinKalGeom_Calo_hh
#include "KinKal/Geometry/Cylinder.hh"
#include "KinKal/Geometry/Disk.hh"
#include "KinKal/Geometry/Annulus.hh"
#include <memory>
namespace mu2e {
  namespace KKGeom {
    class Calo {
      public:
        using CylPtr = std::shared_ptr<KinKal::Cylinder>;
        using DiskPtr = std::shared_ptr<KinKal::Disk>;
        // constructor with geometry parameters from Calorimeter service
        Calo(double z0, double z1, double r0_inner, double r0_outer, double r1_inner, double r1_outer,
             double z0_front, double z0_back, double z1_front, double z1_back);
        // return by reference
        auto const& EMC_Disk_0_Outer() const { return *EMC_Disk_0_Outer_;}
        auto const& EMC_Disk_0_Inner() const { return *EMC_Disk_0_Inner_;}
        auto const& EMC_Disk_1_Inner() const { return *EMC_Disk_1_Inner_;}
        auto const& EMC_Disk_1_Outer() const { return *EMC_Disk_1_Outer_;}

        auto const& EMC_Disk_0_Front() const { return *EMC_Disk_0_Front_;}
        auto const& EMC_Disk_1_Front() const { return *EMC_Disk_1_Front_;}
        auto const& EMC_Disk_0_Back() const  { return *EMC_Disk_0_Back_;}
        auto const& EMC_Disk_1_Back() const  { return *EMC_Disk_1_Back_;}

        // return by ptr
        auto const& EMC_Disk_0_OuterPtr() const { return EMC_Disk_0_Outer_;}
        auto const& EMC_Disk_0_InnerPtr() const { return EMC_Disk_0_Inner_;}
        auto const& EMC_Disk_1_InnerPtr() const { return EMC_Disk_1_Inner_;}
        auto const& EMC_Disk_1_OuterPtr() const { return EMC_Disk_1_Outer_;}
        auto const& EMC_Disk_0_FrontPtr() const { return EMC_Disk_0_Front_;}
        auto const& EMC_Disk_1_FrontPtr() const { return EMC_Disk_1_Front_;}
        auto const& EMC_Disk_0_BackPtr() const  { return EMC_Disk_0_Back_;}
        auto const& EMC_Disk_1_BackPtr() const  { return EMC_Disk_1_Back_;}

        // accessors for local Z positions (relative to tracker center)
        double EMC_Disk_0_Front_Z() const { return z0_front_; }
        double EMC_Disk_0_Back_Z() const { return z0_back_; }
        double EMC_Disk_1_Front_Z() const { return z1_front_; }
        double EMC_Disk_1_Back_Z() const { return z1_back_; }

      private:
        CylPtr EMC_Disk_0_Inner_, EMC_Disk_0_Outer_ , EMC_Disk_1_Inner_, EMC_Disk_1_Outer_;
        DiskPtr EMC_Disk_0_Front_, EMC_Disk_0_Back_, EMC_Disk_1_Front_,  EMC_Disk_1_Back_;
        double z0_front_, z0_back_, z1_front_, z1_back_;  // local Z positions
    };
  }
}

#endif
