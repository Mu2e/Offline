// Original author Jason Bono
//Feb 2018

#include "TrackerConditions/inc/StrawDrift.hh"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include <math.h>
#include <algorithm>
#include "TMath.h"

using namespace std;
namespace mu2e {
  StrawDrift::StrawDrift(std::string filename, float wirevoltage, int phiBins, int dIntegrationBins, float Bz) :
    _phiBins(phiBins)
  {
    //open the file (fixme)
    ifstream myfile(filename);
    if ( !myfile ) {
      //throw cet::exception("FILE");
      cout << "Unable to open particle data file:" << "\n";
    }

    std::vector<double> dataEField;
    std::vector<double> dataVInst;
    std::vector<double> dataDistances; //the set of distances that correspond to a given E-field
    
    //read the file
    string line;
    float a, b;
    while (getline(myfile, line))
    {
      //handle file formating
      istringstream iss(line);
      if ((iss >> a >> b)){
        a = a*100.0; //convert the E-field from KV/cm to V/mm
        b = b*0.01;//convert the speed from cm/us to mm/ns (10/1000 = 0.01)
        //fill the vector of structs
        dataEField.push_back(a);
        dataVInst.push_back(b);
      }
      else {
        //cout << "skipping line"<<"\n";
      }
    }
    
    
    //Use the E:insta-velc tables to build d:insta-veloc tables based on the voltage input
    
    //define the wire and straw radius in mm (remove the hard coding?)
    float wireradius = 12.5/1000.; //12.5 um in mm
    float strawradius = 2.5; //2.5 mm in mm
    
    // calculate the distances that correspond to the efields listed in the table (fix units!!)
    for (size_t i=0; i < dataEField.size(); i++) {
      dataDistances.push_back(wirevoltage/((dataEField[i])*log(strawradius/wireradius))); //in mm
    }
    
    // interpolate to get 50x more points for distances and instantSpeeds
    float fNslices = dIntegrationBins; //for calculations
    float thisDist = 0;
    float nextDist = 0;
    float thisSpeed = 0;
    float nextSpeed = 0;
    float distanceTemp = 0;
    float speedTemp = 0;
    for (size_t i=0; i < (dataDistances.size() - 1); i++) { //cant interpolate beyond the last point (ie avoid i=size+1)
      thisDist = dataDistances[i];
      nextDist = dataDistances[i+1];
      thisSpeed = dataVInst[i];
      nextSpeed = dataVInst[i+1];
      //the linear interpolation
      for (int s=0; s<dIntegrationBins; s++) {
        distanceTemp = ((fNslices - s)*thisDist + (1.0 + s)*nextDist)/(fNslices + 1.0);
        speedTemp = ((fNslices - s)*thisSpeed + (1.0 + s)*nextSpeed)/(fNslices + 1.0);
        this->distances.push_back(distanceTemp);
        this->instantSpeeds.push_back(speedTemp);
      }
    }
    
    //numerically integrate distances and instantSpeeds
    float sliceTime = 0;
    float totalTime = 0;
    float sliceLength = 0;
    float totalLength = 0;
    for (int k = (int) this->distances.size() - 1; k >= 0; k--){ //loop backwards, starting near the wire
      
      if (k == (int) this->distances.size() -1) { //just get the distance to the wire for the closest slice
        sliceLength = this->distances[k];
      }
      else{//get the width of the slices
        sliceLength = this->distances[k] - this->distances[k+1];
      }
      sliceTime = sliceLength/(this->instantSpeeds[k]);
      totalLength += sliceLength;
      totalTime += sliceTime;
      
      this->averageSpeeds.push_back(totalLength/totalTime);
    }
    //Finally, reverse the average speeds vector so that it corresponds with distances; Largest distances first.
    std::reverse(this->averageSpeeds.begin(), this->averageSpeeds.end());
    
    //populate vectors of gamma based on the conditions value of B, and 30 values of phi from 0 to pi/2
    //IN THE GETGAMMA FUNCTION, THE PHI INTERPOLATION WILL BE DONE, AS WELL AS THE 0-2PI ->0-PI/2 REDUNDANCY
    float CC = Bz*log(strawradius/wireradius)/wirevoltage; //multiply this by the average drift velocity to get "C" as defined in doc-5829
    float C = 0;
    float vavg = 0;
    float zetta = 0;
    float dd=0;
    float gammaTemp = 0;
    float vinst = 0;
    int counter = 0;
    //float f_phiBins = _phiBins;
    float phiTemp = 0;
    for (size_t k=0; k < (this->distances.size() - 1); k++) {
      dd = this->distances[k];
      C = CC*vavg*1000000.0;//convert mm/ns to m/s
      vavg = this->averageSpeeds[k];
      zetta = C*dd*0.001;//convert mm to m
      vinst = this->instantSpeeds[k];
      for (int p=0; p < _phiBins; p++) { //this includes both end points (0 and pi/2)
        phiTemp = (float(p)/float(_phiBins - 1.0))*(TMath::Pi()/2.0);
        gammaTemp = (1 + pow(zetta,2)/3)/(1 + pow(zetta*cos(phiTemp),2)/3);
        //fill gamma, effectiveSpeed, time, phi, distance, instantaneousSpeed
        D2Tinfo myD2Tinfo = {gammaTemp, vavg/gammaTemp,dd/(vavg/gammaTemp), phiTemp, dd, vinst};
        this->D2Tinfos.push_back(myD2Tinfo);
        
        counter+=1;
      }
    }
  }//The end of StrawDrift
  
  //look up and return the average speed from vectors
  double StrawDrift::GetAverageSpeed(double dist)
  {
    int lowerIndex = 0;
    float vavg = 0;
    //step through and find distance larger than what is specified
    for (size_t i=0; i < (this->distances.size() - 1); i++) {
      if(dist >= this->distances[i]){
        lowerIndex=i;
        break;
      }
    }
    vavg = this->averageSpeeds[lowerIndex];
    return vavg;
  }
  
  double StrawDrift::GetInstantSpeedFromD(double dist)
  {
    int lowerIndex = 0;
    float vinst = 0;
    //step through and find distance larger than what is specified
    for (size_t i=0; i < (this->distances.size() - 1); i++) {
      if(dist >= this->distances[i]){
        lowerIndex=i;
        break;
      }
    }
    vinst = this->instantSpeeds[lowerIndex];
    return vinst;
  }
  
  double StrawDrift::GetInstantSpeedFromT(double time)
  {
    int lowerIndex = 0;
    float vinst = 0;
    int fullIndex = 0;
    //step through and find time larger than what is specified
    for (size_t i=0; i < (this->distances.size() - 1); i++) {
      fullIndex = i*(_phiBins) + 0; //map from a 2D index to a 1D index (at phi=0)
      if(time >= this->D2Tinfos[fullIndex].time){
        lowerIndex=i;
        break;
      }
    }
    vinst = this->instantSpeeds[lowerIndex];
    return vinst;
  }
  
  double StrawDrift::GetGammaFromD(double distance, double phi)
  {
    float phiSliceWidth = (TMath::Pi()/2.0)/float(_phiBins-1);
    //For the purposes of lorentz corrections, the phi values can be contracted to between 0-90
    float reducedPhi = ConstrainAngle(phi);
    //for interpolation, define a high and a low index
    int upperPhiIndex = ceil(reducedPhi/phiSliceWidth); //rounds the index up to the nearest integer
    int lowerPhiIndex = floor(reducedPhi/phiSliceWidth); //rounds down
    //need the weighting factors
    float lowerPhiWeight = upperPhiIndex - reducedPhi/phiSliceWidth; //a measure of how far the lowerPhiIndex is
    float upperPhiWeight = 1.0 - lowerPhiWeight;
    int fullIndex = 0;
    float upperGamma = 0;
    float lowerGamma = 0;
    for (size_t k=0; k < (this->distances.size() - 1); k++) { //loop through only some of the large vector of structs
      fullIndex = k*(_phiBins)+upperPhiIndex;//mapping from a 2D to a 1D index
      if (distance >= this->D2Tinfos[fullIndex].distance){
        upperGamma = this->D2Tinfos[fullIndex].gamma;//set the gamma associated with the higher index
        fullIndex = k*(_phiBins)+lowerPhiIndex; //mapping from a 2D to a 1D index
        lowerGamma = this->D2Tinfos[fullIndex].gamma;//set the gamma associated with the lower index
        break;
      }
    }
    double Gamma = lowerGamma*lowerPhiWeight + upperGamma*upperPhiWeight;//compute the final gamma
    return Gamma;
  }
  
  double StrawDrift::GetGammaFromT(double time, double phi)
  {
    float phiSliceWidth = (TMath::Pi()/2.0)/float(_phiBins-1);
    //For the purposes of lorentz corrections, the phi values can be contracted to between 0-90
    float reducedPhi = ConstrainAngle(phi);
    //for interpolation, define a high and a low index
    int upperPhiIndex = ceil(reducedPhi/phiSliceWidth); //rounds the index up to the nearest integer
    int lowerPhiIndex = floor(reducedPhi/phiSliceWidth); //rounds down
    //need the weighting factors
    float lowerPhiWeight = upperPhiIndex - reducedPhi/phiSliceWidth; //a measure of how far the lowerPhiIndex is
    float upperPhiWeight = 1.0 - lowerPhiWeight;
    int fullIndex = 0;
    float upperGamma = 0;
    float lowerGamma = 0;
    for (size_t k=0; k < (this->distances.size() - 1); k++) { //loop through only some of the large vector of structs
      fullIndex = k*(_phiBins)+upperPhiIndex;//mapping from a 2D to a 1D index
      if (time >= this->D2Tinfos[fullIndex].time){
        upperGamma = this->D2Tinfos[fullIndex].gamma;//set the gamma associated with the higher index
        fullIndex = k*(_phiBins)+lowerPhiIndex; //mapping from a 2D to a 1D index
        lowerGamma = this->D2Tinfos[fullIndex].gamma;//set the gamma associated with the lower index
        break;
      }
    }
    double Gamma = lowerGamma*lowerPhiWeight + upperGamma*upperPhiWeight;//compute the final gamma
    return Gamma;
  }
  
  //look up and return the lorentz corrected r componenent of the average velocity
  double StrawDrift::GetEffectiveSpeed(double dist, double phi){
    float phiSliceWidth = (TMath::Pi()/2.0)/float(_phiBins-1);
    //For the purposes of lorentz corrections, the phi values can be contracted to between 0-90
    float reducedPhi = fmod(phi,TMath::Pi()/2.0);
    if (reducedPhi < 0){
      reducedPhi += TMath::Pi()/2.0;
    };
    //for interpolation, define a high and a low index
    int upperPhiIndex = ceil(reducedPhi/phiSliceWidth); //rounds the index up to the nearest integer
    int lowerPhiIndex = floor(reducedPhi/phiSliceWidth); //rounds down
    //need the weighting factors
    float lowerPhiWeight = upperPhiIndex - reducedPhi/phiSliceWidth; //a measure of how far the lowerPhiIndex is
    float upperPhiWeight = 1.0 - lowerPhiWeight;
    int fullIndex = 0;
    float upperSpeed = 0;
    float lowerSpeed = 0;
    float effectiveSpeed = 0;
    
    for (size_t k=0; k < (this->distances.size() - 1); k++) { //loop through only some of the large vector of structs
      fullIndex = k*(_phiBins)+upperPhiIndex;
      if (dist >= this->D2Tinfos[fullIndex].distance){
        upperSpeed = this->D2Tinfos[fullIndex].effectiveSpeed; //set the higher speed
        fullIndex = k*(_phiBins)+lowerPhiIndex; // reduce the index by one
        lowerSpeed = this->D2Tinfos[fullIndex].effectiveSpeed; // set the lower speed
        break;
      }
    }
    effectiveSpeed = lowerSpeed*lowerPhiWeight + upperSpeed*upperPhiWeight;
    return effectiveSpeed;
  }
  
  //D2T for sims
  double StrawDrift::D2T(double distance, double phi){
    float phiSliceWidth = (TMath::Pi()/2.0)/float(_phiBins-1);
    //For the purposes of lorentz corrections, the phi values can be contracted to between 0-90
    float reducedPhi = ConstrainAngle(phi);
    //for interpolation, define a high and a low index
    int upperPhiIndex = ceil(reducedPhi/phiSliceWidth); //rounds the index up to the nearest integer
    int lowerPhiIndex = floor(reducedPhi/phiSliceWidth); //rounds down
    //need the weighting factors
    float lowerPhiWeight = upperPhiIndex - reducedPhi/phiSliceWidth; //a measure of how far the lowerPhiIndex is
    float upperPhiWeight = 1.0 - lowerPhiWeight;
    int fullIndex = 0;
    float upperTime = 0;
    float lowerTime = 0;
    float time = 0;
    float gammaTest = 0.;
    for (size_t k=0; k < (this->distances.size() - 1); k++) { //loop through only some of the large vector of structs
      fullIndex = k*(_phiBins)+upperPhiIndex;
      if (distance >= this->D2Tinfos[fullIndex].distance){
        upperTime = this->D2Tinfos[fullIndex].time;//set the higher time
        fullIndex = k*(_phiBins)+lowerPhiIndex; // reduce the index by one
        lowerTime = this->D2Tinfos[fullIndex].time;// set the lower time
        gammaTest = this->D2Tinfos[fullIndex].gamma;//just another test
        break;
      }
    }
    time = lowerTime*lowerPhiWeight + upperTime*upperPhiWeight;//compute the final time
    gammaTest = 1.0*gammaTest;
    return time;
  }
  
  //T2D for reco
  double StrawDrift::T2D(double time, double phi){
    float phiSliceWidth = (TMath::Pi()/2.0)/float(_phiBins-1);
    //For the purposes of lorentz corrections, the phi values can be contracted to between 0-90
    float reducedPhi = fmod(phi,TMath::Pi()/2.0);
    if (reducedPhi < 0){
      reducedPhi += TMath::Pi()/2.0;
    };
    //for interpolation, define a high and a low index
    int upperPhiIndex = ceil(reducedPhi/phiSliceWidth); //rounds the index up to the nearest integer
    int lowerPhiIndex = floor(reducedPhi/phiSliceWidth); //rounds down
    //need the weighting factors
    float lowerPhiWeight = upperPhiIndex - reducedPhi/phiSliceWidth; //a measure of how far the lowerPhiIndex is
    float upperPhiWeight = 1.0 - lowerPhiWeight;
    int fullIndex = 0;
    float upperDist = 0;
    float lowerDist = 0;
    float distance = 0;
    for (size_t k=0; k < (this->distances.size() - 1); k++) { //loop through only some of the large vector of structs
      fullIndex = k*(_phiBins)+upperPhiIndex;
      if (time >= this->D2Tinfos[fullIndex].time){
        upperDist = this->D2Tinfos[fullIndex].distance; //set the higher distance
        fullIndex = k*(_phiBins)+lowerPhiIndex; // reduce the index by one
        lowerDist = this->D2Tinfos[fullIndex].distance;//set the lower distance
        break;
      }
    }
    distance = lowerDist*lowerPhiWeight + upperDist*upperPhiWeight;//compute the final distance
    return distance;
  }
  
  double StrawDrift::ConstrainAngle(double phi){
    if (phi < 0) {
      phi = -1.0*phi;
    }
    phi = fmod(phi,TMath::Pi());
    if (phi > TMath::Pi()/2.0) {
      phi = phi - 2.0*fmod(phi,TMath::Pi()/2.0);
    }
    return phi;
  }
  
}

