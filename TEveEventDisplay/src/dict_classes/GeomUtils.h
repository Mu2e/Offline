#ifndef GeomUtils_h
#define GeomUtils_h
// Math
#include "DataProducts/inc/XYZVec.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "GeometryService/inc/GeomHandle.hh"
//C++
#include <vector>
using namespace CLHEP;

namespace mu2e{
    inline double pointmmTocm(double mm){ return mm/10; };
    inline void hep3vectorTocm(CLHEP::Hep3Vector &vector){vector.set(vector.x()/10, vector.y()/10, vector.z()/10);}
    inline void XYZVecTocm(XYZVec &vector){ vector.SetXYZ(vector.x()/10, vector.y()/10, vector.z()/10);}
   
    inline CLHEP::Hep3Vector GetTrackerCenter(){
      std::string filename("Mu2eG4/geom/geom_common_current.txt");
      SimpleConfig GeomConfig(filename);
      double zCenter  =  GeomConfig.getDouble("mu2e.detectorSystemZ0");
      double xCenter  = -GeomConfig.getDouble("mu2e.solenoidOffset");
      CLHEP::Hep3Vector c(xCenter, 0, zCenter);
      return c;
    }

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

    inline CLHEP::Hep3Vector GetGDMLTrackerOffsetFromMu2e(){ 
      double xCenter = 0;
      double yCenter  =  0;
      double zCenter  = 1288;
      CLHEP::Hep3Vector center(xCenter,yCenter,zCenter);
      return center;
    }

    inline CLHEP::Hep3Vector GetGDMLTrackerCenter() {
      std::string filename("Mu2eG4/geom/mu2eHall.txt");
      SimpleConfig HallConfig(filename);
      double yCenter  = HallConfig.getDouble("yOfFloorSurface.below.mu2eOrigin");
      double zCenter  = 1288; 
      double xCenter  = 0;
      CLHEP::Hep3Vector c(xCenter, yCenter, zCenter);
      return c;
    }

    inline CLHEP::Hep3Vector GetGDMLCaloCenter() {
      std::string filename("Mu2eG4/geom/mu2eHall.txt");
      SimpleConfig HallConfig(filename);
      double yCenter  = HallConfig.getDouble("yOfFloorSurface.below.mu2eOrigin");
      double zCenter  = 1288; 
      double xCenter  = 0;
      CLHEP::Hep3Vector c(xCenter, yCenter, zCenter);
      return c;
    }

    inline CLHEP::Hep3Vector PointToTracker(CLHEP::Hep3Vector point){
      CLHEP::Hep3Vector Mu2eTrackerOrigin = GetTrackerCenter();
      CLHEP::Hep3Vector PointToTracker(point.x() + Mu2eTrackerOrigin.x(), point.y() + Mu2eTrackerOrigin.y(), point.z() + Mu2eTrackerOrigin.z());
      return PointToTracker;
    }

    inline CLHEP::Hep3Vector PointToCalo( CLHEP::Hep3Vector point, int nDisk){
      CLHEP::Hep3Vector Mu2eCaloOrigin = GetCaloCenter(nDisk);
      CLHEP::Hep3Vector PointToCalo(point.x() + Mu2eCaloOrigin.x(), point.y()+Mu2eCaloOrigin.y(), point.z() + Mu2eCaloOrigin.z());
      return  PointToCalo;
    }
    
    inline double TrackerLength(){
      //std::string filename("Mu2eG4/geom/trackerv5.txt");
      //SimpleConfig TrackerConfig(filename);
      //double length  = 2*TrackerConfig.getDouble("tracker.mother.halfLength")/10;
	    return 300.8; //length
    }

    inline double CaloLength(){
      //std::string filename("Mu2eG4/geom/calorimeter_CsI.txt");
      //SimpleConfig CaloConfig(filename);
      //std::vector<double> shift = CaloConfig.getVectorDouble(calorimeter.diskZ0MotherShift); ==700
      //double length  = CaloConfig.getDouble(calorimeter.caloMotherZ1)- CaloConfig.getDouble(calorimeter.caloMotherZ0) - CaloConfig.getDouble(calorimeter.diskCaseZLength); == 1172  
	    return 320;//(shift[1] + length + 1288)/2
    }
}
#endif 
