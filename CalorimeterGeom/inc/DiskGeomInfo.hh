//
// Contains gometry info of disks
//
// Original author B. Echenard
//

#ifndef CalorimeterGeom_DiskGeomInfo_hh
#define CalorimeterGeom_DiskGeomInfo_hh

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"
#include <vector>

namespace mu2e {

    class DiskGeomInfo {

       public:
           DiskGeomInfo() :
             size_                 (CLHEP::Hep3Vector(0,0,0)),
             origin_               (CLHEP::Hep3Vector(0,0,0)),
             originLocal_          (CLHEP::Hep3Vector(0,0,0)),
             rotation_             (CLHEP::HepRotation::IDENTITY),
             inverseRotation_      (CLHEP::HepRotation::IDENTITY),
             originToCrystalOrigin_(CLHEP::Hep3Vector(0,0,0)),
             crystalDirection_     (CLHEP::Hep3Vector(0,0,0)),
             frontFaceCenter_      (CLHEP::Hep3Vector(0,0,0)),
             backFaceCenter_       (CLHEP::Hep3Vector(0,0,0)),
             innerEnvelope_        (0),
             outerEnvelope_        (0),
             crateDeltaZ_          (0)
           {}

           const CLHEP::Hep3Vector&  size()                    const {return size_; }
           const CLHEP::Hep3Vector&  origin()                  const {return origin_;}
           const CLHEP::Hep3Vector&  originLocal()             const {return originLocal_; }
           const CLHEP::Hep3Vector&  originToCrystalOrigin()   const {return originToCrystalOrigin_;}
           const CLHEP::Hep3Vector&  crystalDirection()        const {return crystalDirection_;}
           const CLHEP::HepRotation& rotation()                const {return rotation_;}
           const CLHEP::HepRotation& inverseRotation()         const {return inverseRotation_;}
           const CLHEP::Hep3Vector&  frontFaceCenter()         const {return frontFaceCenter_; }
           const CLHEP::Hep3Vector&  backFaceCenter()          const {return backFaceCenter_; }
           double innerEnvelopeR()                             const {return innerEnvelope_;}
           double outerEnvelopeR()                             const {return outerEnvelope_;}
           double crateDeltaZ()                                const {return crateDeltaZ_;}



           void size(const CLHEP::Hep3Vector& size)                 {size_ = size;}
           void origin(const CLHEP::Hep3Vector& orig)               {origin_ = orig;}
           void originLocal(const CLHEP::Hep3Vector& orig)          {originLocal_ = orig;}
           void originToCrystalOrigin(const CLHEP::Hep3Vector& vec) {originToCrystalOrigin_ = vec;}
           void crystalDirection(const CLHEP::Hep3Vector& vec)      {crystalDirection_ = vec;}
           void rotation(const CLHEP::HepRotation& rot)             {rotation_ = rot; inverseRotation_ = rot.inverse();}
           void frontFaceCenter(const CLHEP::Hep3Vector& pos)       {frontFaceCenter_ = pos;}
           void backFaceCenter(const CLHEP::Hep3Vector& pos)        {backFaceCenter_ = pos;}
           void envelopeRad(double rin, double rout)                {innerEnvelope_ = rin; outerEnvelope_ = rout;}
           void crateDeltaZ(double val)                             {crateDeltaZ_ = val;}


       private:
           CLHEP::Hep3Vector    size_;
           CLHEP::Hep3Vector    origin_;
           CLHEP::Hep3Vector    originLocal_;
           CLHEP::HepRotation   rotation_;
           CLHEP::HepRotation   inverseRotation_;
           CLHEP::Hep3Vector    originToCrystalOrigin_;
           CLHEP::Hep3Vector    crystalDirection_;
           CLHEP::Hep3Vector    frontFaceCenter_;
           CLHEP::Hep3Vector    backFaceCenter_;
           double               innerEnvelope_;
           double               outerEnvelope_;
           double               crateDeltaZ_;
     };
}

#endif
