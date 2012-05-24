#include "ITrackerGeom/inc/CellGeometryHandle_v3.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"

namespace mu2e {

CellGeometryHandle_v3::CellGeometryHandle_v3(ITracker *itr) :
    _itr(itr),
    _nLayer(itr->nSuperLayers())
{
    --_nLayer;
}

CellGeometryHandle_v3::~CellGeometryHandle_v3() {
}

void CellGeometryHandle_v3::SelectCell(int SupLayer, int CelLayer, int Cell, bool isUpstrm) {
        SuperLayer *sl=_itr->getSuperLayer(SupLayer);
        if (SupLayer==_nLayer && CelLayer==1) CelLayer=0;
        if (CelLayer>0)        throw cet::exception("GEOM")<< "The requested cell layer is not allowed in the Square cell geometry"<<std::endl;
        int SelectedCellLayer = 1;
        if (SupLayer==0) SelectedCellLayer++;
        _itl=sl->getLayer(SelectedCellLayer);
        _cell=_itl->getCell(Cell);
        _cell->_tmpMidPoint  = _cell->_senseWire.get()->getMidPoint();
        _cell->_tmpDirection = _cell->_senseWire.get()->getDirection();
        _matrx = _cell->getWire()->get3DTransfrom();
        _invmatrx = _cell->getWire()->get3DInvTransfrom();
        _fSuperLayer=SupLayer;
        _fLayer=CelLayer;
        _fWire=Cell;
        _absRadID = _fSuperLayer;
}

void CellGeometryHandle_v3::SelectCellDet(unsigned long  det) {
        // Return the SuperLayer
        int fSuperLayer=(int)(det*0.0001);

        //Return the Wire
        int fWire=(int)((det)-((fSuperLayer)*10000));

        fSuperLayer--;

        //Call the upper method
        SelectCell(fSuperLayer,0,fWire);
}

void  CellGeometryHandle_v3::SelectCell(int absRadID, int Cell, bool isUpstrm) {
        _absRadID=absRadID;
        SelectCell(absRadID,0,Cell,isUpstrm);
}

unsigned long CellGeometryHandle_v3::computeDet(int SupLayer, int CelLayer, int Cell, bool isUpstrm) {
        unsigned long det = ((SupLayer+1)*10000)+(Cell);
        return det;
}

void  CellGeometryHandle_v3::SelectComp_Cell(int compSupLayer, int compCelLayer, int compICell, bool compIsUpStream) {
        SuperLayer *sl=_itr->getSuperLayer(compSupLayer);
        if (compSupLayer==_nLayer && compCelLayer==1) compCelLayer=0;
        if (compCelLayer>0)        throw cet::exception("GEOM")<< "The requested cell layer is not allowed in the Square cell geometry"<<std::endl;
        int SelectedCellLayer = 1;
        if (compSupLayer==0) SelectedCellLayer++;
        boost::shared_ptr<ITLayer> compItl=sl->getLayer(SelectedCellLayer);
        _Comp_cell = compItl->getCell(compICell);
        _Comp_matrx = _Comp_cell->getWire()->get3DTransfrom();
        _Comp_isUpStream = compIsUpStream;
        _Comp_SuperLayer = compSupLayer;
}

void  CellGeometryHandle_v3::SelectComp_CellDet(unsigned long compDet){
        // Return the SuperLayer
        int compSuperLayer=(int)(compDet*0.0001);

        //Return the Wire
        int compWire=(int)((compDet)-((compSuperLayer)*10000));

        compSuperLayer--;

        SelectComp_Cell(compSuperLayer,0,compWire,false);
}

void  CellGeometryHandle_v3::SelectComp_Cell(int compIAbsRadID, int compICell, bool compIsUpStream) {
        SelectComp_Cell(compIAbsRadID,0,compICell,compIsUpStream);
}

bool  CellGeometryHandle_v3::canIntersectInZ(float &zCorss, float &distWires) const {

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

        maxDistAvail = 0.5*_cell->getWire()->getAlpha() + 0.5*_Comp_cell->getWire()->getAlpha();
        deltaPhi = _cell->_tmpMidPoint.getPhi() - _Comp_cell->_tmpMidPoint.getPhi();

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
                TwoLinePCA pca( _cell->_tmpMidPoint, _cell->_tmpDirection, _Comp_cell->_tmpMidPoint, _Comp_cell->_tmpDirection);
                if (pca.s1()<-_cell->getHalfLength() || pca.s1()>_cell->getHalfLength() ) {
                        return false;
                }
                zCorss = pca.point1().z();//0.5*(pca.point1().z()+pca.point2().z());
                distWires = pca.dca();
        }
        return crossing;
}

bool  CellGeometryHandle_v3::canIntersectInZ(float &zCorss, float &distWires, unsigned long compDet) const {
        // Return the SuperLayer
        int compSuperLayer=(int)(compDet*0.0001);

        //Return the Wire
        int compWire=(int)((compDet)-((compSuperLayer)*10000));

        compSuperLayer--;

        return canIntersectInZ(zCorss,distWires,compSuperLayer,0,compWire);
}

bool  CellGeometryHandle_v3::canIntersectInZ(float &zCorss, float &distWires, int compAbsRadID, int compICell, bool compIsUpstrm) const {
        return canIntersectInZ(zCorss,distWires,compAbsRadID,0,compICell,compIsUpstrm);
}

bool  CellGeometryHandle_v3::canIntersectInZ(float &zCorss, float &distWires, int compSupLayer, int compCelLayer, int compICell, bool compIsUpstrm) const {
        if (_isUpStream!=compIsUpstrm) {
                return false;
        }
        if ( _fSuperLayer == compSupLayer ) {
                return false;
        }
        SuperLayer *sl=_itr->getSuperLayer(compSupLayer);
        if (compSupLayer==_nLayer && compCelLayer==1) compCelLayer=0;
        if (compCelLayer>0)        throw cet::exception("GEOM")<< "The requested cell layer is not allowed in the Square cell geometry"<<std::endl;
        int SelectedCellLayer = 1;
        if (compSupLayer==0) SelectedCellLayer++;
        boost::shared_ptr<ITLayer> compItl=sl->getLayer(SelectedCellLayer);
        boost::shared_ptr<Cell> compCell= compItl->getCell(compICell);

        if (_cell->getWire()->getEpsilon()*compCell->getWire()->getEpsilon()>0.0) {
                return false;
        }

        float deltaPhi, maxDistAvail;
        bool crossing = false;

        maxDistAvail = 0.5*_cell->getWire()->getAlpha() + 0.5*compCell->getWire()->getAlpha();

        deltaPhi = _cell->getMidPoint().getPhi() - compCell->getMidPoint().getPhi();

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
                TwoLinePCA pca(_cell->getMidPoint() , _cell->getDirection(), compCell->getMidPoint(), compCell->getDirection());
                if (pca.s1()<-_cell->getWire()->getDetail()->halfLength() || pca.s1()>_cell->getWire()->getDetail()->halfLength() ) {
                        return false;
                }
                zCorss = pca.point1().z();//0.5*(pca.point1().z()+pca.point2().z());
                distWires = pca.dca();
        }
        return crossing;
}

} // namespace mu2e
