//
// $Id: TrkToCaloExtrapol.cc,v 1.4 2012/09/08 02:24:25 echenard Exp $
// $Author: echenard $
// $Date: 2012/09/08 02:24:25 $
//
// Original author G. Pezzullo
//

// C++ includes
#include <ostream>

// Mu2e includes
#include "TrackCaloMatching/inc/TrkToCaloExtrapol.hh"
#include "CalorimeterGeom/inc/Vane.hh"
#include "CalorimeterGeom/inc/VaneCalorimeter.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "BaBar/BaBar.hh"
#include "TrkBase/HelixParams.hh"
#include "TrkBase/HelixTraj.hh"

using namespace std;

namespace mu2e {




bool TrkToCaloExtrapol::operator == (const TrkToCaloExtrapol & other) const{
        bool res = true;
        res &= (_vaneId == other._vaneId);
        res &= (_pathLengthEntrance == other._pathLengthEntrance);
        res &= (_pathLengthExit == other._pathLengthExit);
        res &= (_trk == other._trk);
        return res;
}
int TrkToCaloExtrapol::vaneId() const{
        return _vaneId;
}

double TrkToCaloExtrapol::time() const{
        return (*_trk.get())->arrivalTime(_pathLengthEntrance);
}

double   TrkToCaloExtrapol::timeErr() const{//FIXME
        return (*_trk.get())->t0().t0Err();
}

double TrkToCaloExtrapol::pathLengthEntrance() const{
        return _pathLengthEntrance;
}

double TrkToCaloExtrapol::t0() const{
        return (*_trk.get())->t0().t0();
}

double TrkToCaloExtrapol::pathLenghtEntranceErr() const{

	GeomHandle<VaneCalorimeter> cg;

        double errTanDip2 = ( (*_trk.get())->helix(_pathLengthEntrance).covariance() )(5, 5);
        double TanDip = (*_trk.get())->helix(_pathLengthEntrance).tanDip();
        HelixTraj trkHel((*_trk.get())->helix(_pathLengthEntrance).params(),
                        (*_trk.get())->helix(_pathLengthEntrance).covariance());
        double z0 = trkHel.z0()
                                        , z = _pathLengthEntrance*trkHel.sinDip() + z0;
        std::cout<<"z0 = "<< z0<<
                        ", z_pathLenght = "<< z<<std::endl;


        Vane const &vane = cg->getVane(_vaneId);
        Hep3Vector Waxes(0.0, 0.0, 1.0);
        Waxes = *(vane.getRotation())*(Waxes);

        Hep3Vector momentumRotUnit = *(vane.getRotation())*TrkToCaloExtrapol::momentum().unit();
        double thetaW = std::atan(-1.0*momentumRotUnit.getZ() / momentumRotUnit.getX() ) ;

        double scaleErrW = 1.0/fabs( cos(thetaW) );

        HepVector momvec(3);
        momvec[0] = Waxes.x();
        momvec[1] = Waxes.y();
        momvec[2] = Waxes.z();

        double tmpErrPos = sqrt(TrkToCaloExtrapol::entrancePositionErr().covMatrix().similarity(momvec));

        double sigmaZ2 = tmpErrPos*scaleErrW;
        std::cout<<"sigmaZ = "<< sigmaZ2<<std::endl;
        sigmaZ2 *= sigmaZ2;

        double res = 0.0;
        double TanDip2 = TanDip*TanDip;
        res += sigmaZ2*(1.0 + TanDip2)/(TanDip2);
        std::cout<<"sigmaS2_z = "<< res << std::endl;

//        double argo = z - z0;
//        argo /= TanDip;
//        argo /= sqrt(1.0 + TanDip);

        double argo = _pathLengthEntrance;
        argo /= ( (1.0 + TanDip2)*TanDip);


        argo = errTanDip2*argo*argo;
        std::cout<<"sigmaS2_tanD = "<< argo << std::endl;
        res += argo;
        return sqrt(res);
}

double    TrkToCaloExtrapol::t0Err() const{
        return (*_trk.get())->t0().t0Err();
}

Hep3Vector TrkToCaloExtrapol::t0Momentum() const{
        return (*_trk.get())->momentum(0.0);
}

BbrVectorErr TrkToCaloExtrapol::t0MomentumErr() const{
        return (*_trk.get())->momentumErr(0.0);
}

double TrkToCaloExtrapol::pathLengthExit() const{
        return _pathLengthExit;
}

HepPoint TrkToCaloExtrapol::entrancePosition() const{
        return (*_trk.get())->position(_pathLengthEntrance);
}

BbrPointErr TrkToCaloExtrapol::entrancePositionErr() const{
        return (*_trk.get())->positionErr(_pathLengthEntrance);
}

HepPoint TrkToCaloExtrapol::exitPosition() const{
        return (*_trk.get())->position(_pathLengthExit);
}

BbrPointErr TrkToCaloExtrapol::exitPositionErr() const{
        return (*_trk.get())->positionErr(_pathLengthExit);
}

Hep3Vector TrkToCaloExtrapol::momentum() const{
        return (*_trk.get())->momentum(_pathLengthEntrance);
}

BbrVectorErr TrkToCaloExtrapol::momentumErr() const{
        return (*_trk.get())->momentumErr(_pathLengthEntrance);
}


// Print the information found in this hit.
void TrkToCaloExtrapol::print( ostream& ost, bool doEndl ) const {

        ost << "TrkToCaloExtrapol :   "
                        << " vane: "          << _vaneId
                        << " time: "          << TrkToCaloExtrapol::time();


        if ( doEndl ){
                ost << endl;
        }

}

} // namespace mu2e
