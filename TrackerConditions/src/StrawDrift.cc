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


using namespace std;
namespace mu2e {
  //StrawDrift::StrawDrift(std::string filename)
    StrawDrift::StrawDrift(std::string filename, float wirevoltage) //JB
  {
    //open the file (fixme)
    ifstream myfile(filename);
    if ( !myfile ) {
        //throw cet::exception("FILE");
        cout << "Unable to open particle data file:" << "\n";
    }
    
    //read the file
    string line;
    float a, b;
    while (getline(myfile, line))
    {
        //handle file formating
        istringstream iss(line);
        if ((iss >> a >> b)){
            //cout << "a = " << a << " b = " <<b<< "\n";
            a = a*100.0; //convert the E-field from KV/cm to V/mm
            b = b*0.01;//convert the speed from cm/us to mm/ns (10/1000 = 0.01)
            //fill the vector of structs
            point mypoint = {a, b};
            this->points.push_back(mypoint);
        }
        else {
            //cout << "skipping line"<<"\n";
        }
    }
    
   // cout << "testing: " << this->points[0].eField<<", "<< this->points[1].eField<<", "<< this->points[2].eField<<", "<< this->points[0].instVel<<", "<< this->points[1].instVel<<", "<< this->points[2].instVel<<"\n";

    
      
      
      //Use the E:insta-velc tables to build d:insta-veloc tables based on the voltage input
      
      //define the wire and straw radius in mm (remove the hard coding?)
      float wireradius = 12.5/1000.; //12.5 um in mm
      float strawradius = 2.5; //2.5 mm in mm
      
      // calculate the distances that correspond to the efields listed in the table (fix units!!)
      for (size_t i=0; i < this->points.size(); i++) {
          this->edistances.push_back(wirevoltage/((points[i].eField)*log(strawradius/wireradius))); //in mm
          cout << "TEST: efield = "<<this->points[i].eField << ", inst-vel = " << this->points[i].instVel << ", edist = " <<this->edistances[i] <<"\n";
          cout << "TEST: wirevoltage = "<<wirevoltage << ", ln = " <<log(strawradius/wireradius)<<"\n";
      }
      

      
      
      
      //cout << "test1" << "\n";
      // interpolate to get 50x more points for distances and instantSpeeds
      const int Nslices = 50;
      float fNslices = Nslices; //for calculations
      float thisDist = 0;
      float nextDist = 0;
      float thisSpeed = 0;
      float nextSpeed = 0;
      float distanceTemp = 0;
      float speedTemp = 0;
      for (size_t i=0; i < (this->points.size() - 1); i++) { //cant interpolate beyond the last point (ie avoid i=size+1)
          thisDist = this->edistances[i];
          nextDist = this->edistances[i+1];
          thisSpeed = this->points[i].instVel;
          nextSpeed = this->points[i+1].instVel;
          //cout <<"testing loop: i="<<i<<"this speed = " <<thisSpeed <<", nextSpeed = " <<nextSpeed <<"\n";
          //the linear interpolation
          //cout << "test2" << "\n";
          for (int s=0; s<Nslices; s++) {
              distanceTemp = ((fNslices - s)*thisDist + (1.0 + s)*nextDist)/(fNslices + 1.0);
              speedTemp = ((fNslices - s)*thisSpeed + (1.0 + s)*nextSpeed)/(fNslices + 1.0);
              this->distances.push_back(distanceTemp);
              this->instantSpeeds.push_back(speedTemp);
          }
      }
      
      //cout << "test3" << "\n";
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
          cout <<"testing loop: k=" <<k<< ", sliceTime ="<<sliceTime<<", sliceLength ="<<sliceLength<<"\n";
          totalLength += sliceLength;
          totalTime += sliceTime;
         
          this->effectiveSpeeds.push_back(totalLength/totalTime);
          //cout << "totalL = " <<totalLength << ", totalTime = " <<totalTime<<", effectiveSpeed = " <<totalLength/totalTime<<"\n";
      }
      //cout << "test4" << "\n";
      //Finally, reverse the effective speeds vector so that it corresponds with distances
      std::reverse(this->effectiveSpeeds.begin(), this->effectiveSpeeds.end());
      
      //cout << "test5" << "\n";

      
      //debugging
      for (int i=0; i < (int) this->distances.size() - 1; i++) {
          cout << "distance = " << this->distances[i] << ", inst V = " <<instantSpeeds[i]<<
          ", effective V = " << this->effectiveSpeeds[i] <<"\n";
      }
      
      
    }
    


    
    //
    //look up and return the effective velocity from vectors
    double StrawDrift::getEffectiveVelocity(double dist)
    {
        int lowerIndex = 0;
        //int upperIndex = 0;
        float vavg = 0;
        //float vlow = 0;
        //float vhigh = 0;
        //step through and first distance larger than what is specified
        for (size_t i=0; i < (this->distances.size() - 1); i++) {
            if(dist >= this->distances[i]){
                lowerIndex=i;
                break;
            }
        }
        vavg = this->effectiveSpeeds[lowerIndex];
        cout << " index="<<lowerIndex<<", dist = " <<dist << ", dist-match = "<< this->distances[lowerIndex] << ", vavg = " << vavg << "\n";
        return vavg;
    }



}

