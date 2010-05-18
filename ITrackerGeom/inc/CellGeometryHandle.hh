#ifndef CELLGEOMETRYHANDLE_HH
#define CELLGEOMETRYHANDLE_HH

#include "CLHEP/Vector/ThreeVector.h"
#include "ITrackerGeom/inc/ITLayer.hh"
#include "ITrackerGeom/inc/Cell.hh"

#include "CLHEP/Geometry/Transform3D.h"
#include "CLHEP/Geometry/Point3D.h"

namespace mu2e {

class CellGeometryHandle {

        friend class ITrackerMaker;

protected:
        CellGeometryHandle();

public:

        ~CellGeometryHandle();

        virtual void  SelectCell(int SupLayer, int CelLayer, int Cell)         {}
        virtual void  SelectWireDet(unsigned long det)                                         {
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
        virtual unsigned long computeDet(int SupLayer, int CelLayer, int Cell) {
                unsigned long det = ((SupLayer+1)*100000)+((CelLayer)*1000)+(Cell);
                return det;
        }
        virtual void  Global2Local(double *global, double *local)                 ;
        virtual void  Local2Global(double *local, double *global)                 ;
        virtual void  WirePosAtEndcap(float *right, float *left)                 ;
        virtual void  WirePosAtZ(float z, float *pos)                                         ;
        virtual void  WirePosAtLength(float length, float *pos)                        ;
        virtual int   GetSuperLayer()                                                                         { return _fSuperLayer; }
        virtual int   GetCelRing()                                                                                 { return _fLayer; }
        virtual int   GetWire()                                                                                 { return _fWire; }
        virtual float GetWireAlfa()                                                                         ;
        virtual float GetWireEpsilon()                                                                         ;
        virtual float GetCellRad()                                                                                 ;
        virtual CLHEP::Hep3Vector GetWireCenter() const                                        ;
        virtual CLHEP::Hep3Vector GetWireDirection() const                                ;

        virtual double DistFromWire(double *global)                                          ;
        virtual double DistFromWireCenter(double *global)                                 ;

        virtual boost::shared_ptr<ITLayer> GetITLayer()                                        { return _itl; }

protected:
        int       _fSuperLayer;
        int       _fLayer;
        int       _fWire;
    boost::shared_ptr<ITLayer> _itl;
    HepGeom::Transform3D _matrx;
    HepGeom::Transform3D _invmatrx;
    boost::shared_ptr<Cell> _cell;
    HepGeom::Point3D<double> tmpGlobal;
    HepGeom::Point3D<double> tmpLocal;
    HepGeom::Point3D<float>  tmpRight;
    HepGeom::Point3D<float>  tmpLeft;
    HepGeom::Point3D<float>  tmpPos;

};

}

#endif /* CELLGEOMETRYHANDLE_HH */
