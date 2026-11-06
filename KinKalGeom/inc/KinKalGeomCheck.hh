#ifndef KinKalGeom_KinKalGeomCheck_hh
#define KinKalGeom_KinKalGeomCheck_hh
//
// Consistency checks of the KinKalGeom passive material: the solid elements it describes (DS shells, passive planes,
// CRV sectors), their pairwise overlaps, coplanar same-material planes, and straight-line path lengths through them,
// used by the KinKalGeomCheck module to compare with the Geant4 (GDML) geometry. All positions in detector coordinates.
//
#include "Offline/KinKalGeom/inc/KinKalGeom.hh"
#include "KinKal/General/Vectors.hh"
#include <random>
#include <string>
#include <vector>

namespace mu2e {
  namespace KKGeomCheck {
    using KinKal::VEC3;
    // A solid element: an oriented box or a tube about wdir; a Gauss-plane pair is merged back into its slab
    struct Element {
      enum Shape { box, tube };
      enum Kind { dsCylinder, passivePlane, crvSector };
      std::string name;
      std::string material;   // KinKal material name; the CRV material for sectors
      Kind kind;
      Shape shape;
      VEC3 center, udir, vdir, wdir; // wdir is the box normal or the tube axis
      double uhalf = 0, vhalf = 0, whalf = 0; // box half lengths
      double rin = 0, rout = 0;               // tube radii (whalf is the tube half length)
      std::vector<size_t> planes;             // passiveMaterialPlanes indices this element came from
      int sector = -1;                        // CRV sector index for a crvSector element
      bool inside(VEC3 const& point) const;
      VEC3 sample(std::mt19937& rng) const;   // uniform in the volume
      double volume() const;
      // length of the segment start + s*dir, s in [smin,smax] (dir a unit vector), inside the element
      double pathLength(VEC3 const& start, VEC3 const& dir, double smin, double smax) const;
      void bounds(VEC3& lo, VEC3& hi) const;  // axis-aligned bounding box
    };
    using Elements = std::vector<Element>;
    // DS material cylinders, passive planes and CRV sectors of a KinKalGeom
    Elements elements(KinKalGeom const& kkg);
    // pairs of elements sharing a volume above 'minVolume' (mm^3), estimated from 'nsample' points per element
    struct Overlap { size_t first, second; double volume; };
    std::vector<Overlap> overlaps(Elements const& elems, unsigned nsample, unsigned seed, double minVolume, bool skipSectorPairs);
    // coplanar CRV sectors with overlapping envelopes: extreme summed layer fractions of lines through the shared area
    struct SectorPair { size_t first, second; double minsum, maxsum; unsigned nlines; };
    std::vector<SectorPair> sectorLayerSums(KinKalGeom const& kkg, double tol, double step);
    // coplanar (within tol) passive rectangles of the same material and thickness whose footprints overlap or touch
    struct Coplanar { size_t first, second; double uoverlap, voverlap; };
    std::vector<Coplanar> coplanarPlanes(KinKalGeom const& kkg, double tol);
    // interval of s in [smin,smax] where start + s*dir is inside a finite cylinder (axis through 'center'); false if none
    bool cylinderInterval(VEC3 const& center, VEC3 const& axis, double radius, double halflen,
        VEC3 const& start, VEC3 const& dir, double smin, double smax, double& s0, double& s1);
  }
}
#endif
