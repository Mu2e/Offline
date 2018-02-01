#include "TrackerConditions/inc/StrawDrift.hh"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>


using namespace std;
namespace mu2e {
  StrawDrift::StrawDrift(std::string filename)
  {
    //open the file
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
            cout << "a = " << a << " b = " <<b<< "\n";
            
            //fill the vector of structs
            point mypoint = {a, b};
            this->points.push_back(mypoint);
        }
        else {
            cout << "skipping line"<<"\n";
        }
    }
    
    cout << "testing: " << this->points[0].eField<<", "<< this->points[1].eField<<", "<< this->points[2].eField<<", "<< this->points[0].instVel<<", "<< this->points[1].instVel<<", "<< this->points[2].instVel<<"\n";
  }
}

