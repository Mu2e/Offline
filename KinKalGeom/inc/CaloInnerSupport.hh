//
//  Inner support of a calorimeter disk crystal case: three solid aluminium rings (front, middle and back) and a thin
//  carbon fiber ring between them and the crystals, see caloBuildCase in Mu2eG4/src/constructDiskCalorimeter.cc.
//  A track reaches this material either through the upstream face of the case, just behind the front panel, or from
//  the disk hole, through the inner cylinder of the case. It is described by the (r,z) cross sections of the rings,
//  and the material crossed is the straight-line path through them from the entry point until the track reaches the
//  crystals, returns to the hole or leaves the case. Once a track has reached the crystals it showers, so only the
//  entry into the case is considered.
//  original author: Sophie Middleton (2026)
//
#ifndef KinKalGeom_CaloInnerSupport_hh
#define KinKalGeom_CaloInnerSupport_hh
#include "Offline/DataProducts/inc/SurfaceId.hh"
#include "KinKal/General/Vectors.hh"
#include "KinKal/Geometry/Annulus.hh"
#include "KinKal/Geometry/Cylinder.hh"
#include "KinKal/Geometry/Plane.hh"
#include <array>
#include <memory>
#include <string>
#include <vector>
namespace mu2e {
  namespace KKGeom {
    class CaloInnerSupport {
      public:
        using AnnPtr = std::shared_ptr<KinKal::Annulus>;
        using CylPtr = std::shared_ptr<KinKal::Cylinder>;
        using PlanePtr = std::shared_ptr<KinKal::Plane>;
        enum Material { aluminium=0, carbon, nmaterials };
        // the ways a track can cross the support material, used as the surface Id index:
        // aluminium entered through the case front or from the hole, the carbon fiber ring entered from its inner
        // radius, or end-on through the case front, and aluminium entered from the carbon fiber ring inner radius
        enum Crossing { alFront=0, alHole, cfRing, cfFront, alRing, ncrossings };
        // ring cross section in (r,z), z in detector coordinates
        struct Ring {
          double rin_, rout_, zmin_, zmax_;
          Material mat_;
        };
        // the case spans [zfront,zback] between radius rin and the crystals at rcrystal; the carbon fiber ring
        // starts at rcf
        CaloInnerSupport(SurfaceIdEnum::enum_type sid, KinKal::VEC3 const& center, double rin, double rcf, double rcrystal,
            double zfront, double zback, std::array<std::string,nmaterials> const& materials);
        void addRing(Ring const& ring) { rings_.push_back(ring); }
        SurfaceId surfaceId(Crossing xing) const { return SurfaceId(sid_,int(xing)); }
        static Material material(Crossing xing) { return xing == alFront || xing == alHole || xing == alRing ? aluminium : carbon; }
        auto const& materialName(Material mat) const { return materials_[mat]; }
        auto const& front() const { return front_; } // case front between the hole and the crystals
        auto const& frontPlane() const { return frontplane_; } // unbounded plane of the case front, to find where a track reaches it
        auto const& inner() const { return inner_; } // case inner cylinder
        auto const& cfInner() const { return cfinner_; } // carbon fiber ring inner cylinder
        auto const& rings() const { return rings_; }
        double rin() const { return inner_->radius(); }
        double rcf() const { return cfinner_->radius(); }
        double rcrystal() const { return rcrystal_; }
        double transverseRadius(KinKal::VEC3 const& pos) const { return (pos - center_).Rho(); }
        // straight-line path length (mm) through each material from an entry point on the case surface in the given
        // direction, until the track reaches the crystals, goes back into the hole or leaves the case
        std::array<double,nmaterials> pathLengths(KinKal::VEC3 const& pos, KinKal::VEC3 const& dir) const;
      private:
        SurfaceIdEnum::enum_type sid_;
        KinKal::VEC3 center_; // disk axis position at the case front
        double rcrystal_;
        PlanePtr frontplane_;
        AnnPtr front_;
        CylPtr inner_;
        CylPtr cfinner_;
        std::array<std::string,nmaterials> materials_; // KinKal material names
        std::vector<Ring> rings_;
    };
  }
}
#endif
