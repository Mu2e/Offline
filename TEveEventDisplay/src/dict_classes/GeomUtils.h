#ifndef GeomUtils_h
#define GeomUtils_h
// Math
#include "Offline/DataProducts/inc/GenVector.hh"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
//C++
#include <vector>
using namespace CLHEP;

namespace mu2e{
    inline double pointmmTocm(double mm){ return mm/10; };
    inline void hep3vectorTocm(CLHEP::Hep3Vector &vector){vector.set(vector.x()/10, vector.y()/10, vector.z()/10);}
    inline void XYZVectorFTocm(XYZVectorF &vector){ vector.SetXYZ(vector.x()/10, vector.y()/10, vector.z()/10);}
   
   inline CLHEP::Hep3Vector GetCaloCenter(int nDisk){
      std::string calfilename("Mu2eG4/geom/calorimeter_CsI.txt");
      SimpleConfig CalConfig(calfilename);
      double zCenter = 0;
      if(nDisk==0) zCenter = CalConfig.getDouble("calorimeter.caloMotherZ0") + 100;
      if(nDisk==1) zCenter = CalConfig.getDouble("calorimeter.caloMotherZ1") - 600;
      std::string geomfilename("Mu2eG4/geom/geom_common_current.txt");
      SimpleConfig GeomConfig(geomfilename);
      double xCenter  = -GeomConfig.getDouble("mu2e.solenoidOffset");
      CLHEP::Hep3Vector c(xCenter, 0, zCenter);
      return c;
    }


    inline CLHEP::Hep3Vector PointToCaloCM( CLHEP::Hep3Vector point, int nDisk){
      CLHEP::Hep3Vector Mu2eCaloOrigin = GetCaloCenter(nDisk);
      CLHEP::Hep3Vector PointToCalo(10*(point.x() + Mu2eCaloOrigin.x()), 10*(point.y()+Mu2eCaloOrigin.y()), 10*(point.z() + Mu2eCaloOrigin.z()));
      return  PointToCalo;
    }
    
    inline CLHEP::Hep3Vector PointToCalo( CLHEP::Hep3Vector point, int nDisk){
      CLHEP::Hep3Vector Mu2eCaloOrigin = GetCaloCenter(nDisk);
      CLHEP::Hep3Vector PointToCalo((point.x() + Mu2eCaloOrigin.x()), (point.y()+Mu2eCaloOrigin.y()), (point.z() + Mu2eCaloOrigin.z()));
      return  PointToCalo;
    }

    inline double TrackerLength(){
      GeomHandle<Tracker> trkr;
      TubsParams envelope(trkr->g4Tracker()->getInnerTrackerEnvelopeParams());
      double dz{(envelope.zHalfLength())};
      return (dz*2);
    }
    
    inline CLHEP::Hep3Vector NewCenter(){
      GeomHandle<Tracker> trkr;
      GeomHandle<DetectorSystem> det;
      CLHEP::Hep3Vector origin = trkr->origin();
      CLHEP::Hep3Vector InMu2e = det->toMu2e(origin);
      return InMu2e;
    }

    inline double CaloLength(){
      //Not used anymore....
      return 320;
    }
}
#endif 
