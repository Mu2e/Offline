#include "ITrackerGeom/inc/CellGeometryHandle.hh"

#include "cetlib/pow.h"

using cet::square;

namespace mu2e {

void CellGeometryHandle::Global2Local(double *global, double *local)
{
        tmpGlobal.set(global[0],global[1],global[2]);
        tmpLocal=_invmatrx*tmpGlobal;
        local[0]=tmpLocal.x();
        local[1]=tmpLocal.y();
        local[2]=tmpLocal.z();

//        CLHEP::Hep3Vector tr  = _matrx.getTranslation();
//        double mt0  = global[0]-tr[0];
//        double mt1  = global[1]-tr[1];
//        double mt2  = global[2]-tr[2];
//        CLHEP::HepRotation rot = _matrx.getRotation();
//        local[0] = mt0*rot[0][0] + mt1*rot[1][0] + mt2*rot[2][0];
//        local[1] = mt0*rot[0][1] + mt1*rot[1][1] + mt2*rot[2][1];
//        local[2] = mt0*rot[0][2] + mt1*rot[1][2] + mt2*rot[2][2];

}

void CellGeometryHandle::Local2Global(double *local, double *global)
{
        tmpLocal.set(local[0],local[1],local[2]);
        tmpGlobal=_matrx*tmpLocal;
        global[0]=tmpGlobal.x();
        global[1]=tmpGlobal.y();
        global[2]=tmpGlobal.z();

//        int i;
//        CLHEP::Hep3Vector tr  = _matrx.getTranslation();
//        CLHEP::HepRotation rot = _matrx.getRotation();
//        for (i=0; i<3; i++) {
//                global[i] = tr[i]
//                            + local[0]*rot[i][0]
//                            + local[1]*rot[i][1]
//                            + local[2]*rot[i][2];
//        }

}

void CellGeometryHandle::WirePosAtEndcap(float *right, float *left)
{
        tmpRight.set(0.0,0.0,_cell->getWire()->getDetail()->halfLength());
        tmpLeft.set(0.0,0.0,-tmpRight.z());
        tmpRight.transform(_matrx);
        tmpLeft.transform(_matrx);
        right[0]=tmpRight.x();
        right[1]=tmpRight.y();
        right[2]=tmpRight.z();
        left[0]=tmpLeft.x();
        left[1]=tmpLeft.y();
        left[2]=tmpLeft.z();
}

void CellGeometryHandle::WirePosAtZ(float z, float *pos)
{
        WirePosAtLength(z/cos(_cell->getWire()->getEpsilon()), pos);
}

void CellGeometryHandle::WirePosAtLength(float length, float *pos)
{
        tmpPos.set(0.0,0.0,length);
        tmpPos.transform(_matrx);
        pos[0]=tmpPos.x();
        pos[1]=tmpPos.y();
        pos[2]=tmpPos.z();
}

float CellGeometryHandle::GetWireAlfa()
{
        return (float) _cell->getWire()->getAlpha();
}

float CellGeometryHandle::GetWireEpsilon()
{
        return (float) _cell->getWire()->getEpsilon();
}

float CellGeometryHandle::GetCellRad()
{
        return (float) _cell->getDetail()->CirumscribedRadius();
}

CLHEP::Hep3Vector CellGeometryHandle::GetWireCenter() const {
        return _cell->getMidPoint();
}

CLHEP::Hep3Vector CellGeometryHandle::GetWireDirection() const {
        return _cell->getDirection();
}

double CellGeometryHandle::DistFromWire(double *global)
{
        return (DistFromWireCenter(global)-_cell->getDetail()->wireRadius());
}

double CellGeometryHandle::DistFromWireCenter(double *global)
{
        tmpGlobal.set(global[0],global[1],global[2]);
        tmpLocal=_invmatrx*tmpGlobal;
        return sqrt(square(tmpLocal.x())+square(tmpLocal.y()));
}

} // namespace mu2e
