// ITracker wire description
//
// $Id: Wire.cc,v 1.7 2012/09/25 10:08:28 tassiell Exp $
// $Author: tassiell $
// $Date: 2012/09/25 10:08:28 $
//
// Original author G. Tassielli
//

#include "ITrackerGeom/inc/Wire.hh"
#include "CLHEP/Geometry/Vector3D.h"

#ifndef __CINT__

namespace mu2e {

Wire::Wire():
  _id(WireId()),
  _wireType(undefined),
  _c(CLHEP::Hep3Vector(0.,0.,0.)),
  _w(CLHEP::Hep3Vector(0.,0.,1.)),
  _pos(0x0),
  _invpos(HepGeom::Transform3D()),
  _epsilon(0.0),
  _alpha(0.0)
{
}


Wire::Wire( WireId id,
                boost::shared_ptr<WireDetail> detail,
                HepGeom::Transform3D *pos,
                double epsilon,
                double alpha,
                Wtype wireType
            ):
  _id(id),
  _wireType(wireType),
  _pos(pos),
  _epsilon(epsilon),
  _alpha(alpha)
{
  _detail=detail;
  _c.set(_pos->dx(), _pos->dy(), _pos->dz());
  _invpos=_pos->inverse();
  _w.set(_pos->xz(),_pos->yz(),_pos->zz());

}

Wire::~Wire (){
//        try {
//                delete *_detail; *_detail=NULL;
//                delete _pos;
//        } catch (cet::exception e) {
//            throw cet::exception("GEOM")
//                << "Error during deleting wire data \n";
//        }

  delete _pos;
}

} // namespace mu2e

#endif
