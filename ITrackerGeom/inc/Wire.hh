#ifndef ITrackerGeom_Wire_hh
#define ITrackerGeom_Wire_hh

#include <deque>
#include <vector>

#include <boost/shared_ptr.hpp>

#include "ITrackerGeom/inc/WireDetail.hh"
#include "ITrackerGeom/inc/WireId.hh"

#include "CLHEP/Geometry/Transform3D.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

class Wire{

        friend class ITracker;
        friend class ITrackerMaker;


public:

        // A free function, returning void, that takes a const Wire& as an argument.
        typedef void (*WireFunction)( const Wire& s);

        enum Wtype {undefined=-1, field, sense};

        Wire();

        // Constructor using wire 3d transformation.
        Wire( WireId id,
                        boost::shared_ptr<WireDetail> detail,
                        HepGeom::Transform3D *pos,
                        double epsilon,
                        double alpha,
                        Wtype wireType=undefined
        );

        ~Wire ();

        WireId Id() const { return _id;}

        boost::shared_ptr<WireDetail> getDetail() const { return _detail;}

        HepGeom::Transform3D get3DTransfrom() const {return *_pos;}

        HepGeom::Transform3D get3DInvTransfrom() const {return _invpos;}

        CLHEP::Hep3Vector getMidPoint() const {return _c;}

        CLHEP::Hep3Vector getDirection() const { return _w;}

        double getEpsilon() const { return _epsilon;}
        double getAlpha() const { return _alpha;}

        Wtype getWireType() const { return _wireType; }

protected:

        // Identifier
        WireId _id;

        Wtype _wireType;

        // Mid-point of the wire.
        CLHEP::Hep3Vector _c;

        // Detailed description of a wire.
        boost::shared_ptr<WireDetail> _detail;

        // Unit vector along the wire direction.
        // Need to add unit vectors along local u and v also.
        // Use Euler angle convention from G4.
        CLHEP::Hep3Vector _w;

        const HepGeom::Transform3D *_pos;
        HepGeom::Transform3D _invpos;
        const double _epsilon;
        const double _alpha;

};

inline std::ostream& operator<<(std::ostream& ost, const Wire& w ){
        ost <<w.Id()<<" type "<<w.getWireType()<<" radius "<< w.getDetail()->outerRadius() << " length "<< w.getDetail()->length() <<std::endl;
        ost<<"epsilon "<< w.getEpsilon()<<" pos matrix: "<<std::endl;
        ost<<w.get3DTransfrom().xx()<<" "<<w.get3DTransfrom().xy()<<" "<<w.get3DTransfrom().xz()<<" "<<w.get3DTransfrom().dx()<<std::endl;
        ost<<w.get3DTransfrom().yx()<<" "<<w.get3DTransfrom().yy()<<" "<<w.get3DTransfrom().yz()<<" "<<w.get3DTransfrom().dy()<<std::endl;
        ost<<w.get3DTransfrom().zx()<<" "<<w.get3DTransfrom().zy()<<" "<<w.get3DTransfrom().zz()<<" "<<w.get3DTransfrom().dz()<<std::endl;
        return ost;
}

}  //namespace mu2e

#endif /* ITrackerGeom_Wire_hh */
