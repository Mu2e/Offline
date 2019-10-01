//
// Hold information about one Straw.
//
// Original author Rob Kutschke
//

#include <sstream>
#include "TrackerGeom/inc/Straw.hh"
#include "TrackerGeom/inc/Tracker.hh"

#ifndef __CINT__

using std::vector;
using CLHEP::Hep3Vector;

namespace mu2e {

  Straw::Straw():
    _id(StrawId()),
    _c(CLHEP::Hep3Vector(0.,0.,0.)),
    _w(CLHEP::Hep3Vector(0.,0.,1.))
  {
  }

  Straw::Straw( const StrawId& id,
                CLHEP::Hep3Vector const& c,
                double wtx,
                double wty
                ):
    _id(id),
    _c(c)
  {
    _w = CLHEP::Hep3Vector(wtx,wty,1.).unit();
  }

  Straw::Straw( const StrawId& id,
                CLHEP::Hep3Vector const& c,
                CLHEP::Hep3Vector const& w
                ):
    _id(id),
    _c(c),
    _w(w.unit())
  {
  }

  void Straw::fillPointers ( const Tracker* tracker ) const{
    _tracker = tracker;
  }

  std::string Straw::name( std::string const& base ) const{
    std::ostringstream os;

    os << base
       << _id.getPlane() << "_"
       << _id.getPanel() << "_"
       << _id.getLayer()  << "_"
       << _id.getStraw();
    return os.str();
  }

    // outer radius
   double Straw::getRadius() const { return _tracker->strawOuterRadius(); }

   double Straw::innerRadius() const { return _tracker->strawInnerRadius(); }

   double  Straw::getThickness() const { return _tracker->strawWallThickness(); }

    // Half length
   double Straw::halfLength() const { 
               return _tracker->getStrawHalfLength(_id.straw()); }

    // active length is a little smaller
   double Straw::activeHalfLength() const { 
             return _tracker->getStrawActiveHalfLength(_id.straw()); }



} // namespace mu2e
#endif
