#include "TRandom3.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
// generate parameters to misalign the tracker to include an overall twist (dalphaz/dz), skew (dalphay/dz), and random mis-positions and angles
// David Brown (LBNL)
//

using namespace std;
void MisalignTracker(double twist, double skew, double squeeze, double sigalpha, double sigpos, const char* outfile) {
  unsigned NPlanes(36);
  double planezgap(57.0); // FIXME!!
  
  TRandom3 myrand(3499803);

  filebuf fb;
  fb.open (outfile,ios::out);
  ostream os(&fb);
  // write header
  os << "# Tracker Misalignments from MisalignTracker.C, parameters: "
  << " twist = "<< twist
  << " skew = "<< skew
  << " sigalpha = "<< sigalpha
  << " sigpos = "<< sigpos
  << endl;
  for(unsigned iplane=0;iplane < NPlanes; ++iplane){
    double planez = iplane*planezgap;
    double dalphax = myrand.Gaus(0.0,sigalpha);
    double dalphay = myrand.Gaus(skew,sigalpha);
    double dalphaz = myrand.Gaus(twist*planez,sigalpha);
    os << "Plane " << iplane 
      << " dalpha : [ "
      << dalphax << " , "
      << dalphay << " , "
      << dalphaz << " ] "
      << " dpos : [ "
      << myrand.Gaus(0.0,sigpos) << " , "
      << myrand.Gaus(0.0,sigpos) << " , "
      << myrand.Gaus(0.0,sigpos) << " ] "
      << endl;
  }

}

