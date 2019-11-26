#include "TRandom3.h"
#include "TMath.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
// generate parameters to misalign the tracker to include an overall twist (dphi/zdz), skew in x and y, and random mis-positions and angles
// David Brown (LBNL)
//

using namespace std;
void MisalignTracker(double twist, double skew, double squeeze, double sigalpha, double sigpos, const char* outfile) {
  unsigned NPlanes(36);
  double planezgap(83.0); // approximate number
  double planer(700.0); // approximate number
  double trackerlen = NPlanes*planezgap; // approximate

  TRandom3 myrand(3499803);
  double rtwist = myrand.Uniform(-twist,twist);
  double xskew = myrand.Uniform(-skew,skew);
  double yskew = myrand.Uniform(-skew,skew);
  double rsqueeze = myrand.Uniform(-squeeze,squeeze);
  double dtwist = twist/(TMath::Pi()*2.0*planer*trackerlen);

  filebuf fb;
  fb.open (outfile,ios::out);
  ostream os(&fb);
  // write header
  os << "# Tracker Misalignments from MisalignTracker.C, parameters: "
  << " twist = "<< twist << " mm" 
  << " skew = "<< skew << " mm" 
  << " squeeze = "<< squeeze << " mm" 
  << " rtwist = "<< rtwist << " mm" 
  << " dtwist = " << dtwist << " rad/mm"
  << " xskew = "<< xskew << " mm" 
  << " yskew = "<< yskew << " mm" 
  << " rsqueeze = "<< rsqueeze << " mm" 
  << " sigalpha = "<< sigalpha << " rad" 
  << " sigpos = "<< sigpos << " mm"
  << endl;
  for(unsigned ibuf=0;ibuf<3;ibuf++)
    os << "#" << endl;

  // global track: dummy for now
  os << " TABLE TrkAlignTracker " << endl;
  os << "0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0" << endl;
  for(unsigned ibuf=0;ibuf<3;ibuf++)
    os << "#" << endl;
 
  os << "TABLE TrkAlignPlane" << endl;
  for(unsigned iplane=0;iplane < NPlanes; ++iplane){
    double planez = (iplane-float(NPlanes)/2.0)*planezgap;
    double dalphax = myrand.Gaus(0.0,sigalpha);
    double dalphay = myrand.Gaus(0.0,sigalpha);
    double dalphaz = myrand.Gaus(dtwist*planez,sigalpha);
    double dx = myrand.Gaus(planez*xskew/trackerlen,sigpos);
    double dy = myrand.Gaus(planez*yskew/trackerlen,sigpos);
    double dz = myrand.Gaus(planez*rsqueeze/trackerlen,sigpos);

    os << iplane << ", " 
      << dx << ", "
      << dy << ", "
      << dz << ", "
      << dalphax << ", "
      << dalphay << ", "
      << dalphaz << endl;
  }

   for(unsigned ibuf=0;ibuf<3;ibuf++)
    os << "#" << endl;
 // dummy panel alignment for now
  os << "TABLE TrkAlignPanel" << endl;
  for(unsigned ipanel=0;ipanel<6*NPlanes; ++ipanel){
    os << ipanel << ", "
    << "0.0, 0.0, 0.0, 0.0, 0.0, 0.0" << endl;
  }
}
