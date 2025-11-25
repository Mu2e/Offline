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
        auto const& EMC_Disk_0_SurfIn() const { return *EMC_Disk_0_SurfIn_; }
        auto const& EMC_Disk_0_SurfOut() const { return *EMC_Disk_0_SurfOut_; }
        auto const& EMC_Disk_1_SurfIn() const { return *EMC_Disk_1_SurfIn_; }
        auto const& EMC_Disk_1_SurfOut() const  { return *EMC_Disk_1_SurfOut_;  }

        auto const& EMC_Disk_0_EdgeIn() const { return *EMC_Disk_0_EdgeIn_; }
        auto const& EMC_Disk_0_EdgeOut() const { return *EMC_Disk_0_EdgeOut_; }
        auto const& EMC_Disk_1_EdgeIn() const { return *EMC_Disk_1_EdgeIn_; }
        auto const& EMC_Disk_1_EdgeOut() const  { return *EMC_Disk_1_EdgeOut_;  }

        auto const& EMC_0_FrontIn() const { return *EMC_0_FrontIn_; }
        auto const& EMC_0_FrontOut() const { return *EMC_0_FrontOut_; }
        auto const& EMC_1_FrontIn() const { return *EMC_1_FrontIn_; }
        auto const& EMC_1_FrontOut() const  { return *EMC_1_FrontOut_;  }

        auto const& EMC_2_FrontIn() const { return *EMC_2_FrontIn_; }
        auto const& EMC_2_FrontOut() const { return *EMC_2_FrontOut_; }
        auto const& EMC_3_FrontIn() const { return *EMC_3_FrontIn_; }
        auto const& EMC_3_FrontOut() const  { return *EMC_3_FrontOut_;  }

// return by ptr
        auto const& EMC_Disk_0_SurfInPtr() const { return EMC_Disk_0_SurfIn_;}
        auto const& EMC_Disk_0_SurfOutPtr() const { return EMC_Disk_0_SurfOut_;}
        auto const& EMC_Disk_1_SurfInPtr() const { return EMC_Disk_1_SurfIn_;}
        auto const& EMC_Disk_1_SurfOutPtr() const  { return EMC_Disk_1_SurfOut_;}

        auto const& EMC_Disk_0_EdgeInPtr() const { return EMC_Disk_0_EdgeIn_;}
        auto const& EMC_Disk_0_EdgeOutPtr() const { return EMC_Disk_0_EdgeOut_;}
        auto const& EMC_Disk_1_EdgeInPtr() const { return EMC_Disk_1_EdgeIn_;}
        auto const& EMC_Disk_1_EdgeOutPtr() const  { return EMC_Disk_1_EdgeOut_;}

        auto const& EMC_0_FrontInPtr() const { return EMC_0_FrontIn_;}
        auto const& EMC_0_FrontOutPtr() const { return EMC_0_FrontOut_;}
        auto const& EMC_1_FrontInPtr() const { return EMC_1_FrontIn_;}
        auto const& EMC_1_FrontOutPtr() const  { return EMC_1_FrontOut_;}

        auto const& EMC_2_FrontInPtr() const { return EMC_2_FrontIn_;}
        auto const& EMC_2_FrontOutPtr() const { return EMC_2_FrontOut_;}
        auto const& EMC_3_FrontInPtr() const { return EMC_3_FrontIn_;}
        auto const& EMC_3_FrontOutPtr() const  { return EMC_3_FrontOut_;}

      private:
        CylPtr EMC_Disk_0_SurfIn_, EMC_Disk_0_SurfOut_, EMC_Disk_1_SurfIn_, EMC_Disk_1_SurfOut_, EMC_Disk_0_EdgeIn_, EMC_Disk_0_EdgeOut_,EMC_Disk_1_EdgeIn_, EMC_Disk_1_EdgeOut_; // active volume boundary in XY amd YZ

        DiskPtr EMC_0_FrontIn_, EMC_0_FrontOut_, EMC_1_FrontIn_, EMC_1_FrontOut_, EMC_2_FrontIn_, EMC_2_FrontOut_, EMC_3_FrontIn_, EMC_3_FrontOut_;

    };
  }
}

#endif
