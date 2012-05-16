#include "ITrackerGeom/inc/CellGeometryHandle_v2.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"

namespace mu2e {

CellGeometryHandle_v2::CellGeometryHandle_v2(ITracker *itr) {
        _itr=itr;
}

CellGeometryHandle_v2::~CellGeometryHandle_v2() {
}

void CellGeometryHandle_v2::SelectCell() {
        SuperLayer *sl=_itr->getSuperLayer(_fSuperLayer);
        int SelectedCellLayer = _fLayer*2;
        if (_fLayer==0) SelectedCellLayer+=2;
        if (_fSuperLayer==0) SelectedCellLayer++;
        _itl=sl->getLayer(SelectedCellLayer);
        _cell=_itl->getCell(_fWire);
        _cell->_tmpMidPoint  = _cell->_senseWire.get()->getMidPoint();
        _cell->_tmpDirection = _cell->_senseWire.get()->getDirection();
 //        _cell = _itr->getSuperLayer(_fSuperLayer)->getLayer(CelLayer)->getCell(Cell);
        _matrx = _cell->getWire()->get3DTransfrom();
        _invmatrx = _cell->getWire()->get3DInvTransfrom();
}

void CellGeometryHandle_v2::SelectCell(int SupLayer, int CelLayer, int Cell, bool isUpstrm) {
        _fSuperLayer=SupLayer;
        _fLayer=CelLayer;
        _fWire=Cell;
        _absRadID = SupLayer*_itr->nRing()+CelLayer;
        SelectCell();
}

void  CellGeometryHandle_v2::SelectCell(int absRadID, int Cell, bool isUpstrm) {
        _absRadID=absRadID;
        _fSuperLayer=absRadID/_itr->nRing();
        _fLayer=absRadID-_itr->nRing()*_fSuperLayer;
        _fWire=Cell;
        SelectCell();
}

int CellGeometryHandle_v2::computeAbsRadID(int SupLayer, int CelLayer, bool isUpstrm) {
        int absRadID = SupLayer*_itr->nRing()+CelLayer;
        return absRadID;
}

bool  CellGeometryHandle_v2::canIntersectInZ(float &zCorss, float &distWires, int compAbsRadID, int compICell, bool compIsUpstrm) const {
        int compSupLayer=compAbsRadID/_itr->nRing();
        int compCelLayer=compAbsRadID-_itr->nRing()*compSupLayer;
        return canIntersectInZ(zCorss,distWires,compSupLayer,compCelLayer,compICell,compIsUpstrm);
}

bool  CellGeometryHandle_v2::canIntersectInZ(float &zCorss, float &distWires, int compSupLayer, int compCelLayer, int compICell, bool compIsUpstrm) const {
        if (_isUpStream!=compIsUpstrm) {
                return false;
        }
        if ( _fSuperLayer == compSupLayer ) {
                return false;
        }

        SuperLayer *sl=_itr->getSuperLayer(compSupLayer);
        int SelectedCellLayer = compCelLayer*2;
        if (compCelLayer==0) SelectedCellLayer+=2;
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
                zCorss = pca.point1().z();
                distWires = pca.dca();
        }
        return crossing;
}

} // namespace mu2e
