// Stack of the planes containing sensor modules
//
// Evan Schiewe, 2013

#ifndef EXTMONFNALPLANESTACK_HH
#define EXTMONFNALPLANESTACK_HH

#include <vector>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALPlane.hh"
#include "canvas/Persistency/Common/Wrapper.h"

namespace mu2e {

  namespace ExtMonFNAL { class ExtMonMaker; }
  namespace ExtMonFNAL { class ExtMon; }

  class ExtMonFNALPlaneStack {
  public:

    unsigned nplanes() const { return m_plane_zoffset.size(); }
    unsigned nModulesPerPlane() const { return planes_[0].module_zoffset().size(); }
    unsigned nmodules() const { return (this->nplanes() * this->nModulesPerPlane()); }

    const std::vector<ExtMonFNALPlane>& planes() const { return planes_; }
    const int size() const { return planes_.size(); }

    const std::vector<double>& plane_zoffset() const { return m_plane_zoffset; }
    const std::vector<double>& plane_xoffset() const { return m_plane_xoffset; }
    const std::vector<double>& plane_yoffset() const { return m_plane_yoffset; }

    // offset of plane center wrt the ref point
    CLHEP::Hep3Vector   planeOffsetInStack(unsigned iplane) const;

    // Positioning of the stack: reference point and rotation.
    CLHEP::Hep3Vector refPointInMu2e() const { return m_stackRefPointInMu2e; }
    CLHEP::HepRotation const& rotationInMu2e() const { return m_stackRotationInMu2e; }

    // Coordinate conversion to/from the Mu2e frame
    // The ExtMonFNAL frame is defined in the following way:
    //
    // - The (0,0,0) point is the reference point of the stack,
    //
    // - The z_stack axis is perpendicular to the plane planes
    // - The x_stack axis in in the horizontal plane
    // - The y_stack axis forms a right-handed (x_em, y_em, z_em) frame
    //
    // The rotation w.r.t. Mu2e is "small", that is, the projection of
    // z_stack on z_mu2e is positive and similar for x, y.

    CLHEP::Hep3Vector mu2eToStack_position(const CLHEP::Hep3Vector& mu2epos) const;
    CLHEP::Hep3Vector mu2eToStack_momentum(const CLHEP::Hep3Vector& mu2emom) const;

    CLHEP::Hep3Vector stackToMu2e_position(const CLHEP::Hep3Vector& pos) const;
    CLHEP::Hep3Vector stackToMu2e_momentum(const CLHEP::Hep3Vector& mom) const;


    //----------------------------------------------------------------
    // "global" extmon plane number is obtained by adding the offset to
    // this stack's plane number
    unsigned planeNumberOffset() const { return planeNumberOffset_; }

    //----------------------------------------------------------------
  private:
    ExtMonFNALPlaneStack();
    friend class ExtMonFNAL::ExtMon;
    friend class ExtMonFNAL::ExtMonMaker;

    // For persistency
    template<class T> friend class art::Wrapper;

    std::vector<ExtMonFNALPlane> planes_;


    ExtMonFNALPlane plane_;

    CLHEP::HepRotation m_stackRotationInMu2e;
    CLHEP::Hep3Vector m_stackRefPointInMu2e;

    // Mother volume (G4Polyhedra)
    std::vector<double> m_motherTransverseHalfSize;
    double m_motherStartZ;
    double m_motherEndZ;

    // Plane center positions
    std::vector<double> m_plane_zoffset;
    std::vector<double> m_plane_xoffset;
    std::vector<double> m_plane_yoffset;

    // data for coordinate system transformations: inverse of stack rotation
    CLHEP::HepRotation m_coordinateRotationInMu2e;

    unsigned planeNumberOffset_;
  };

  //================================================================

} // namespace mu2e

#endif/*EXTMONFNALPLANESTACK_HH*/
