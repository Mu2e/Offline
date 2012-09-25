// interface to manage the geometries of the ITracker cells
//
// $Id: CellGeometryHandle_v3_DBL.cc,v 1.3 2012/09/25 10:08:28 tassiell Exp $
// $Author: tassiell $
// $Date: 2012/09/25 10:08:28 $
//
// Original author G. Tassielli
//

#include "ITrackerGeom/inc/CellGeometryHandle_v3_DBL.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"

namespace mu2e {

CellGeometryHandle_v3_DBL::CellGeometryHandle_v3_DBL(ITracker *itr):CellGeometryHandle_v3(itr){
        _nSLayer=itr->nSuperLayers();
}


CellGeometryHandle_v3_DBL::~CellGeometryHandle_v3_DBL() {
}

inline void  CellGeometryHandle_v3_DBL::AdjustMatrix() {
        if (_isDownStream) {
                _DnStrmDeadWireLngt=_itr->zZonesLimits()[1]/cos(_cell->getWire()->getEpsilon());
                _tmpCellHalfLength=_cell->getHalfLength()-_DnStrmDeadWireLngt;
                _tmpCellHalfLength*=0.5;
                WirePosAtLength(_tmpCellHalfLength+_DnStrmDeadWireLngt,wCntPos);
                _tmpMidPoint.set(wCntPos[0],wCntPos[1],wCntPos[2]);
                _tmpDirection.set(_cell->getDirection().x(),_cell->getDirection().y(),_cell->getDirection().z());
        }
        else if (_isUpStream) {
                _UpStrmDeadWireLngt=_itr->zZonesLimits()[0]/cos(_cell->getWire()->getEpsilon());
                _tmpCellHalfLength=_cell->getHalfLength()+_UpStrmDeadWireLngt;
                _tmpCellHalfLength*=0.5;
                WirePosAtLength(-_tmpCellHalfLength+_UpStrmDeadWireLngt,wCntPos);
                _tmpMidPoint.set(wCntPos[0],wCntPos[1],wCntPos[2]);
                _tmpDirection.set(-_cell->getDirection().x(),-_cell->getDirection().y(),-_cell->getDirection().z());
        }
}

void CellGeometryHandle_v3_DBL::SelectCell(int SupLayer, int CelLayer, int Cell, bool isUpstrm) {
        if (SupLayer>=_nSLayer || isUpstrm){
                _isUpStream=true;
                _isDownStream=false;
                SupLayer=SupLayer%_nSLayer;
        }
        else {
                _isUpStream=false;
                _isDownStream=true;
        }
        CellGeometryHandle_v3::SelectCell(SupLayer,CelLayer,Cell);
        AdjustMatrix();
}

void CellGeometryHandle_v3_DBL::SelectCellDet(unsigned long  det) {
        // Return the SuperLayer
        int fSuperLayer=(int)(det*0.0001);

        //Return the Wire
        int fWire=(int)((det)-((fSuperLayer)*10000));

        fSuperLayer--;
        if (fSuperLayer>=_nSLayer){
                _isUpStream=true;
                _isDownStream=false;
                fSuperLayer=fSuperLayer%_nSLayer;
        }
        else {
                _isUpStream=false;
                _isDownStream=true;
        }

        //Call the upper method
        CellGeometryHandle_v3::SelectCell(fSuperLayer,0,fWire);
        AdjustMatrix();
}

void  CellGeometryHandle_v3_DBL::SelectCell(int absRadID, int Cell, bool isUpstrm) {
        SelectCell(absRadID,0,Cell,isUpstrm);
}

unsigned long CellGeometryHandle_v3_DBL::computeDet(int SupLayer, int CelLayer, int Cell, bool isUpstrm) {
        if (SupLayer<_nSLayer && isUpstrm) {
                SupLayer+=_nSLayer;
        }
        return CellGeometryHandle_v3::computeDet(SupLayer,CelLayer,Cell);
}

int CellGeometryHandle_v3_DBL::computeAbsRadID(int SupLayer, int CelLayer, bool isUpstrm) {
        if (SupLayer>=_nSLayer || isUpstrm){
                SupLayer=SupLayer%_nSLayer;
        }
        int absRadID = SupLayer;
        return absRadID;
}

const CLHEP::Hep3Vector& CellGeometryHandle_v3_DBL::GetCellCenter() const {
        return _tmpMidPoint;
}

const CLHEP::Hep3Vector& CellGeometryHandle_v3_DBL::GetCellDirection() const {
        return _tmpDirection;
}

double CellGeometryHandle_v3_DBL::GetCellHalfLength() const {
        return _tmpCellHalfLength;
}

void CellGeometryHandle_v3_DBL::Comp_WirePosAtLength(float length, float *pos)
{
        HepGeom::Point3D<float> Comp_tmpPos;
        Comp_tmpPos.set(0.0,0.0,length);
        Comp_tmpPos.transform(_Comp_matrx);
        pos[0]=Comp_tmpPos.x();
        pos[1]=Comp_tmpPos.y();
        pos[2]=Comp_tmpPos.z();
}

inline void  CellGeometryHandle_v3_DBL::Comp_AdjustMatrix() {
        if (_Comp_isDownStream) {
                _Comp_DnStrmDeadWireLngt=_itr->zZonesLimits()[1]/cos(_Comp_cell->getWire()->getEpsilon());
                _Comp_tmpCellHalfLength=_Comp_cell->getHalfLength()-_Comp_DnStrmDeadWireLngt;
                _Comp_tmpCellHalfLength*=0.5;
                Comp_WirePosAtLength(_Comp_tmpCellHalfLength+_Comp_DnStrmDeadWireLngt,Comp_wCntPos);
                _Comp_tmpMidPoint.set(Comp_wCntPos[0],Comp_wCntPos[1],Comp_wCntPos[2]);
                _Comp_tmpDirection.set(_Comp_cell->getDirection().x(),_Comp_cell->getDirection().y(),_Comp_cell->getDirection().z());
        }
        else if (_Comp_isUpStream) {
                _Comp_UpStrmDeadWireLngt=_itr->zZonesLimits()[0]/cos(_Comp_cell->getWire()->getEpsilon());
                _Comp_tmpCellHalfLength=_Comp_cell->getHalfLength()+_Comp_UpStrmDeadWireLngt;
                _Comp_tmpCellHalfLength*=0.5;
                Comp_WirePosAtLength(-_Comp_tmpCellHalfLength+_Comp_UpStrmDeadWireLngt,Comp_wCntPos);
                _Comp_tmpMidPoint.set(Comp_wCntPos[0],Comp_wCntPos[1],Comp_wCntPos[2]);
                _Comp_tmpDirection.set(-_Comp_cell->getDirection().x(),-_Comp_cell->getDirection().y(),-_Comp_cell->getDirection().z());
        }
}

void  CellGeometryHandle_v3_DBL::SelectComp_Cell(int compSupLayer, int compCelLayer, int compICell, bool compIsUpStream) {
        if (compSupLayer>=_nSLayer || compIsUpStream){
                _Comp_isUpStream=true;
                _Comp_isDownStream=false;
                compSupLayer=compSupLayer%_nSLayer;
         }
         else {
                 _Comp_isUpStream=false;
                 _Comp_isDownStream=true;
         }
        CellGeometryHandle_v3::SelectComp_Cell( compSupLayer, compCelLayer, compICell, compIsUpStream);
        Comp_AdjustMatrix();
}

void  CellGeometryHandle_v3_DBL::SelectComp_CellDet(unsigned long compDet){
        // Return the SuperLayer
        int compSuperLayer=(int)(compDet*0.0001);

        //Return the Wire
        int compWire=(int)((compDet)-((compSuperLayer)*10000));

        compSuperLayer--;
        bool compIsUpStream = false;
        if (compSuperLayer>=_nSLayer){
                compIsUpStream=true;
                compSuperLayer=compSuperLayer%_nSLayer;
        }
        else {
                compIsUpStream=false;
        }

        SelectComp_Cell(compSuperLayer,0,compWire,compIsUpStream);
}

void  CellGeometryHandle_v3_DBL::SelectComp_Cell(int compIAbsRadID, int compICell, bool compIsUpStream) {
        SelectComp_Cell(compIAbsRadID,0,compICell,compIsUpStream);
}

bool  CellGeometryHandle_v3_DBL::canIntersectInZ(float &zCorss, float &distWires) const {

        if (_isUpStream!=_Comp_isUpStream) {
                return false;
        }
        if ( _fSuperLayer == _Comp_SuperLayer ) {
                return false;
        }

        if (_cell->getWire()->getEpsilon()*_Comp_cell->getWire()->getEpsilon()>0.0) {
                return false;
        }

        float deltaPhi, maxDistAvail;
        bool crossing = false;

        maxDistAvail = 0.5*_tmpCellHalfLength/_cell->getHalfLength()*_cell->getWire()->getAlpha() +
                        0.5*_Comp_tmpCellHalfLength/_Comp_cell->getHalfLength()*_Comp_cell->getWire()->getAlpha();
        deltaPhi = _tmpMidPoint.getPhi() - _Comp_tmpMidPoint.getPhi();

        if ( deltaPhi==0.0 ) {
                crossing = true;
        } else {
                deltaPhi = fabs (deltaPhi);
                if (deltaPhi>CLHEP::pi) deltaPhi = CLHEP::twopi - deltaPhi;
                if (deltaPhi<maxDistAvail) {
                        crossing = true;
                }
        }
        if ( crossing ) {
                TwoLinePCA pca( _tmpMidPoint, _tmpDirection, _Comp_tmpMidPoint, _Comp_tmpDirection);
                if (pca.s1()<-_tmpCellHalfLength || pca.s1()>_tmpCellHalfLength ) {
                        return false;
                }
                zCorss = pca.point1().z();//0.5*(pca.point1().z()+pca.point2().z());
                distWires = pca.dca();
        }
        return crossing;
}

} // namespace mu2e
