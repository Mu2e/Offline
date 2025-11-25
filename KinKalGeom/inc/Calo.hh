//
//  Define the nominal calorimeter boundary and reference surfaces, used to extrapolate and sample KinKal track fits, and to build
//  the passive materials in the fit
//  original author: Sophie Middleton (2025)
//
#ifndef KinKalGeom_Calorimeter_hh
#define KinKalGeom_Calorimeter_hh
#include "KinKal/Geometry/Cylinder.hh"
#include "KinKal/Geometry/Disk.hh"
#include "KinKal/Geometry/Annulus.hh"
#include <memory>
namespace mu2e {
  namespace KinKalGeom {
    class Calorimeter {
      public:
        using CylPtr = std::shared_ptr<KinKal::Cylinder>;
        using DiskPtr = std::shared_ptr<KinKal::Disk>;
        // default constructor with nominal geometry
        Calorimeter();
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

      private:
        CylPtr EMC_Disk_0_Outer, EMC_Disk_0_Inner, EMC_Disk_1_Inner, EMC_Disk_1_Outer;
        DiskPtr EMC_Disk_0_Front, EMC_Disk_1_Front, EMC_Disk_0_Back, EMC_Disk_1_Back;
    };
  }
}

#endif
