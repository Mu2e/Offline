#ifndef CalorimeterGeom_DiskInfo_hh
#define CalorimeterGeom_DiskInfo_hh
//
// Contains geometry info of a disk.
//
// Design: the disk POSE (origin_ + rotation_) and the body-fixed LOCAL offsets
// are the only source of truth. Every global-frame point (face centers, crystal
// direction) and inverseRotation_ are DERIVED and cached, rebuilt in recompute_()
// whenever the pose or a local offset changes. This makes rigid moves correct by
// construction: moveDisk only needs to update the pose.
//
// NOTE for the maker: set the pose (origin + rotation) BEFORE the frontFaceCenter/
// backFaceCenter/crystalDirection setters, since those back-solve the stored local
// offset from the global value using the current pose.
//
// Original author B. Echenard
//

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"
#include <iostream>

namespace mu2e {

  class DiskInfo {

    public:
      DiskInfo();

      // ---- frame conversions (single source of the transform) ----
      CLHEP::Hep3Vector toGlobal(const CLHEP::Hep3Vector& local)  const;
      CLHEP::Hep3Vector toLocal (const CLHEP::Hep3Vector& global) const;
      CLHEP::Hep3Vector toLocalFF (const CLHEP::Hep3Vector& g) const;
      CLHEP::Hep3Vector toGlobalFF(const CLHEP::Hep3Vector& l) const;

      // ---- getters (derived points return the cached value) ----
      const CLHEP::Hep3Vector&  size()                  const {return size_;}
      const CLHEP::Hep3Vector&  origin()                const {return origin_;}
      const CLHEP::Hep3Vector&  originLocal()           const {return originLocal_;}
      const CLHEP::Hep3Vector&  originToCrystalOrigin() const {return originToCrystalOrigin_;}
      const CLHEP::Hep3Vector&  crystalDirection()      const {return crystalDirection_;}
      const CLHEP::Hep3Vector&  frontFaceCenter()       const {return frontFaceCenter_;}
      const CLHEP::Hep3Vector&  backFaceCenter()        const {return backFaceCenter_;}
      const CLHEP::HepRotation& rotation()              const {return rotation_;}
      const CLHEP::HepRotation& inverseRotation()       const {return inverseRotation_;}
      double crystalZLength()                           const {return crystalZlength_;}
      double innerEnvelopeR()                           const {return innerEnvelope_;}
      double outerEnvelopeR()                           const {return outerEnvelope_;}
      double FEBZOffset()                               const {return FEBZOffset_;}
      double FEBZLength()                               const {return FEBZLength_;}

      // ---- pose: the only mutable geometric state ----
      void origin  (const CLHEP::Hep3Vector& o)  {origin_ = o;   recompute_();}
      void rotation(const CLHEP::HepRotation& r) {rotation_ = r; recompute_();}
      void setPose (const CLHEP::Hep3Vector& o, const CLHEP::HepRotation& r)
                   {origin_ = o; rotation_ = r; recompute_();}

      // ---- body-fixed / invariant state ----
      void originLocal          (const CLHEP::Hep3Vector& o) {originLocal_ = o;}
      void originToCrystalOrigin(const CLHEP::Hep3Vector& v) {originToCrystalOrigin_ = v;}   // local, invariant under moves

      // ---- global-valued setters: back-solve the local offset, then rebuild ----
      void frontFaceCenter (const CLHEP::Hep3Vector& g) {ffLocal_  = toLocal(g);  recompute_();}
      void backFaceCenter  (const CLHEP::Hep3Vector& g) {bfLocal_  = toLocal(g);  recompute_();}
      void crystalDirection(const CLHEP::Hep3Vector& d) {dirLocal_ = rotation_*d; recompute_();}   // direction: rotate only

      // ---- frame-invariant scalars / packed sizes ----
      void size          (const CLHEP::Hep3Vector& s) {size_ = s;}
      void crystalZlength(double v)                   {crystalZlength_ = v;}
      void envelopeRad   (double rin, double rout)    {innerEnvelope_ = rin; outerEnvelope_ = rout;}
      void FEBZOffset    (double v)                   {FEBZOffset_ = v;}
      void FEBZLength    (double v)                   {FEBZLength_ = v;}

      void print(std::ostream& os = std::cout) const;


    private:
      void recompute_();

      // pose (source of truth)
      CLHEP::Hep3Vector  origin_{0,0,0};
      CLHEP::Hep3Vector  originLocal_{0,0,0};
      CLHEP::HepRotation rotation_{CLHEP::HepRotation::IDENTITY};

      // body-fixed offsets (source of truth, invariant under rigid moves)
      CLHEP::Hep3Vector  size_{0,0,0};
      CLHEP::Hep3Vector  originToCrystalOrigin_{0,0,0};
      CLHEP::Hep3Vector  ffLocal_{0,0,0};
      CLHEP::Hep3Vector  bfLocal_{0,0,0};
      CLHEP::Hep3Vector  dirLocal_{0,0,0};

      // derived caches (rebuilt only in recompute_)
      CLHEP::HepRotation inverseRotation_{CLHEP::HepRotation::IDENTITY};
      CLHEP::Hep3Vector  frontFaceCenter_{0,0,0};
      CLHEP::Hep3Vector  backFaceCenter_{0,0,0};
      CLHEP::Hep3Vector  crystalDirection_{0,0,0};

      double crystalZlength_{0};
      double innerEnvelope_{0};
      double outerEnvelope_{0};
      double FEBZOffset_{0};
      double FEBZLength_{0};
  };
}

#endif
