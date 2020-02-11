#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include "TMath.h"

#include "TrackerConditions/inc/StrawDriftMaker.hh"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "GeometryService/inc/GeomHandle.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "BFieldGeom/inc/BFieldManager.hh"
#include "BTrk/BField/BField.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "CLHEP/Matrix/Vector.h"


using namespace std;

namespace mu2e {

  StrawDrift::ptr_t StrawDriftMaker::fromFcl() {

    double wirevoltage = _config.wireVoltage();

    std::vector<double> dataEField = _config.kVcm();
    std::vector<double> dataVInst = _config.cmus();
    if ( dataEField.size()==0 || dataEField.size()!=dataVInst.size() ) {
      throw cet::exception("STRAW_DRIFT_BADMODEL")
	<< "input drift model don't make sense, sizes:"<< dataEField.size() 
	<< " " << dataVInst.size()  << "\n";
    }
    //the set of distances that correspond to a given E-field
    std::vector<double> dataDistances; 
    //convert the E-field from KV/cm to V/mm
    for(auto& x : dataEField) x *= 100.0;
    //convert the speed from cm/us to mm/ns (10/1000 = 0.01)
    for(auto& x : dataVInst) x *= 0.01;

    // Use the E:insta-velc tables to build d:insta-veloc tables 
    // based on the voltage input
    
    //define the wire and straw radius in mm
    GeomHandle<Tracker> tracker;
    //12.5 um in mm
    double wireradius = tracker->wireRadius(); 
    //2.5 mm in mm
    double strawradius = tracker->strawOuterRadius(); 
    
    // calculate the distances that correspond to the efields 
    // listed in the table (fix units!!)
    double logRadius = log(strawradius/wireradius);
    for (size_t i=0; i < dataEField.size(); i++) {
      dataDistances.push_back(wirevoltage/((dataEField[i])*logRadius)); //in mm
    }

    // quantities to be saved:
    std::vector<float> distances; // 
    std::vector<float> instantSpeeds; // the instantaneous "nominal" speed
    std::vector<float> averageSpeeds; // the average "nominal" speed
    std::vector<StrawDrift::D2Tinfo> D2Tinfos; // 2-D array in distance and phi

    // interpolate to get driftIntegrationBins more points 
    // for distances and instantSpeeds
    float fNslices = _config.driftIntegrationBins(); //for calculations
    float thisDist = 0;
    float nextDist = 0;
    float thisSpeed = 0;
    float nextSpeed = 0;
    float distanceTemp = 0;
    float speedTemp = 0;
    //cant interpolate beyond the last point (ie avoid i=size+1)
    for (size_t i=0; i < (dataDistances.size() - 1); i++) { 
      thisDist = dataDistances[i];
      nextDist = dataDistances[i+1];
      thisSpeed = dataVInst[i];
      nextSpeed = dataVInst[i+1];
      //the linear interpolation
      for (int s=0; s<fNslices; s++) {
        distanceTemp = ((fNslices - s)*thisDist + (1.0 + s)*nextDist)/(fNslices + 1.0);
        speedTemp = ((fNslices - s)*thisSpeed + (1.0 + s)*nextSpeed)/(fNslices + 1.0);
        distances.push_back(distanceTemp);
        instantSpeeds.push_back(speedTemp);
      }
    }
    
    //numerically integrate distances and instantSpeeds
    float sliceTime = 0;
    float totalTime = 0;
    float sliceLength = 0;
    float totalLength = 0;
    //loop backwards, starting near the wire
    for (int k = (int) distances.size() - 1; k >= 0; k--){   
      if (k == (int) distances.size() -1) { 
	//just get the distance to the wire for the closest slice
        sliceLength = distances[k];
      } else {
	//get the width of the slices
        sliceLength = distances[k] - distances[k+1];
      }
      sliceTime = sliceLength/instantSpeeds[k];
      totalLength += sliceLength;
      totalTime += sliceTime;
      
      averageSpeeds.push_back(totalLength/totalTime);
    }
    //Finally, reverse the average speeds vector so that it 
    // corresponds with distances; Largest distances first.
    std::reverse(averageSpeeds.begin(), averageSpeeds.end());

    size_t phiBins = _config.phiBins();

    // populate vectors of gamma based on the conditions value of B, 
    // and phiBins values of phi from 0 to pi/2
    //IN THE GETGAMMA FUNCTION, THE PHI INTERPOLATION WILL BE DONE, AS WELL AS THE 0-2PI ->0-PI/2 REDUNDANCY
    GeomHandle<BFieldManager> bfmgr;
    GeomHandle<DetectorSystem> det;
    CLHEP::Hep3Vector vpoint_mu2e = det->toMu2e(CLHEP::Hep3Vector(0.0,0.0,0.0));
    float Bz = bfmgr->getBField(vpoint_mu2e).z();

    float CC = Bz*logRadius/wirevoltage; 
    //multiply this by the average drift velocity to get "C" as defined in doc-5829
    float C = 0;
    float vavg = 0;
    float zetta = 0;
    float dd=0;
    float gammaTemp = 0;
    float vinst = 0;
    float phiTemp = 0;
    for (size_t k=0; k < (distances.size() - 1); k++) {
      dd = distances[k];
      C = CC*vavg*1000000.0;//convert mm/ns to m/s
      vavg = averageSpeeds[k];
      zetta = C*dd*0.001;//convert mm to m
      vinst = instantSpeeds[k];
      for (size_t p=0; p < phiBins; p++) { //this includes both end points (0 and pi/2)
        phiTemp = (float(p)/float(phiBins - 1.0))*(TMath::Pi()/2.0);
        gammaTemp = (1 + pow(zetta,2)/3)/(1 + pow(zetta*cos(phiTemp),2)/3);
        //fill gamma, effectiveSpeed, time, phi, distance, instantaneousSpeed
	StrawDrift::D2Tinfo myD2Tinfo = {gammaTemp, vavg/gammaTemp,
			     dd/(vavg/gammaTemp), phiTemp, dd, vinst};
	D2Tinfos.push_back(myD2Tinfo);        
      }
    }

    auto ptr = std::make_shared<StrawDrift>(D2Tinfos,distances,
					    instantSpeeds,averageSpeeds,phiBins);
    
    return ptr;

  } // fromFcl
  

  StrawDrift::ptr_t StrawDriftMaker::fromDb() {
    return fromFcl();
  }

}

