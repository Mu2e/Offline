#ifndef ITrackerGeom_CellGeometryHandle_hh
#define ITrackerGeom_CellGeometryHandle_hh

#include "CLHEP/Vector/ThreeVector.h"
#include "ITrackerGeom/inc/Cell.hh"
#include "ITrackerGeom/inc/ITLayer.hh"

#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Geometry/Transform3D.h"

namespace mu2e {

class CellGeometryHandle {

        friend class ITrackerMaker;

protected:
        CellGeometryHandle() {
                _isDownStream=true;
                _isUpStream=false;
        }

public:

        virtual ~CellGeometryHandle() {}

        virtual void  SelectCell(int SupLayer, int CelLayer, int Cell, bool isUpstrm=false) {
                throw cet::exception("GEOM")<< "SelectCell Method not implemented for the interface class CellGeometryHandle, please use one of the real implementation"<<std::endl;
        }
        virtual void  SelectCellDet(unsigned long det) {
                // Return the SuperLayer
                int fSuperLayer=(int)(det*0.00001);

                //Return the Layer
                int fLayer=(int)((det)-((fSuperLayer)*100000));

                fLayer=(int)(fLayer*0.001);

                //Return the Wire
                int fWire=(int)(((det)-((fSuperLayer)*100000))-fLayer*1000);

                fSuperLayer--;

                //Call the upper method
                SelectCell(fSuperLayer,fLayer,fWire);
        }
        virtual void  SelectCell(int absRadID, int Cell, bool isUpstrm=false) {
                _absRadID=absRadID;
                SelectCell(absRadID,0,Cell,isUpstrm);
        }
        virtual unsigned long computeDet(int SupLayer, int CelLayer, int Cell, bool isUpstrm=false) {
                unsigned long det = ((SupLayer+1)*100000)+((CelLayer)*1000)+(Cell);
                return det;
        }
        virtual int computeAbsRadID(int SupLayer, int CelLayer, bool isUpstrm=false) {
                int absRadID = SupLayer;
                return absRadID;
        }
        virtual void  Global2Local(double *global, double *local) ;
        virtual void  Local2Global(double *local, double *global) ;
        virtual void  WirePosAtEndcap(float *right, float *left)  ;
        virtual void  WirePosAtZ(float z, float *pos)             ;
        virtual void  WirePosAtLength(float length, float *pos)   ;
        virtual void  WirePosAtEndcap(CLHEP::Hep3Vector &right, CLHEP::Hep3Vector &left)  ;
        virtual void  WirePosAtZ(float z, CLHEP::Hep3Vector &pos)             ;
        virtual void  WirePosAtLength(float length, CLHEP::Hep3Vector &pos)   ;
        virtual int   GetSuperLayer()                             { return _fSuperLayer; }
        virtual int   GetCelRing()                                { return _fLayer; }
        virtual int   GetWire()                                   { return _fWire; }
        virtual float GetWireAlfa()                               ;
        virtual float GetWireEpsilon()                            ;
        virtual float GetCellRad()                                ;
        virtual int   GetCellAbsRadID()                           { return _absRadID; }
        virtual const CLHEP::Hep3Vector& GetWireCenter() const    ;
        virtual const CLHEP::Hep3Vector& GetWireDirection() const ;
        virtual const CLHEP::Hep3Vector& GetCellCenter() const    ;
        virtual const CLHEP::Hep3Vector& GetCellDirection() const ;
        virtual double GetCellHalfLength() const                  ;

        virtual double DistFromWire(double *global)               ;
        virtual double DistFromWireCenter(double *global)         ;
        virtual double DistFromWire(CLHEP::Hep3Vector &global)            ;
        virtual double DistFromWireCenter(CLHEP::Hep3Vector &global)      ;
        virtual double DistFromWire(CLHEP::Hep3Vector const &global)      ;
        virtual double DistFromWireCenter(CLHEP::Hep3Vector const &global);

        virtual boost::shared_ptr<ITLayer> GetITLayer()           { return _itl; }
        virtual boost::shared_ptr<Cell> GetITCell()               { return _cell; }
        bool isUpStream()   { return _isUpStream; }
        bool isDownStream() { return _isDownStream; }

protected:
        int       _fSuperLayer;
        int       _fLayer;
        int       _fWire;
        int       _absRadID;
        boost::shared_ptr<ITLayer> _itl;
        HepGeom::Transform3D _matrx;
        HepGeom::Transform3D _invmatrx;
        boost::shared_ptr<Cell> _cell;
        HepGeom::Point3D<double> tmpGlobal;
        HepGeom::Point3D<double> tmpLocal;
        HepGeom::Point3D<float>  tmpRight;
        HepGeom::Point3D<float>  tmpLeft;
        HepGeom::Point3D<float>  tmpPos;

        bool _isUpStream;
        bool _isDownStream;

};

}

#endif /* ITrackerGeom_CellGeometryHandle_hh */
