//
// generate parameters to misalign the tracker to include an overall twist (dphi/zdz), skew in x and y, and random mis-positions and angles
// David Brown (LBNL)
//
#include "TRandom3.h"
#include "TMath.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>

using namespace std;

  // panel azimuth
double PanelPhi(int plane, int panel){
  double pphi0[6] { 106, 195, 225,315, 345,75};
  double pphi1[6] { 255, 165, 135, 45, 15,285};
  int mplane = plane%4;
  double radfac =M_PI/180;
  if(mplane==0 || mplane ==3) {
    return radfac*pphi0[panel];
  } else {
    return radfac*pphi1[panel];
  }
}


void MisalignTracker(int seed, double weakfactor, double trackerangle, double trackerpos,
    double planeangle, double planepos,
    double panelangle, double panelpos,
    double wiredV, double wiredW,
    const char* outfile="MisalignTracker.txt") {
  unsigned NPlanes(36);
  unsigned NPanels(6);
  unsigned NStraws(96);
  double planezgap(83.0); // approximate number
  double planer(700.0); // approximate number
  double trackerlen = NPlanes*planezgap; // approximate

  TRandom3 myrand(seed);
  // define weak modes from tracker parameters
  double weakangle = trackerangle*weakfactor;
  double weakpos = trackerpos*weakfactor;
  double xskew = myrand.Uniform(-weakangle,weakangle);
  double yskew = myrand.Uniform(-weakangle,weakangle);
  double twist = myrand.Uniform(-weakangle,weakangle);
  double zsqueeze = myrand.Uniform(-weakpos,weakpos);
  double rsqueeze = myrand.Uniform(-weakpos,weakpos);

  filebuf fb;
  fb.open (outfile,ios::out);
  ostream os(&fb);
  // write header
  os << "# Tracker Misalignments from MisalignTracker.C, parameters: " << endl
    << "# weakfactor " << weakfactor <<  endl
    << "# trackerangle " << trackerangle << " rad" << endl
    << "# trackerpos " << trackerpos <<" mm" << endl
    << "# planeangle = "<< planeangle << " rad" << endl
    << "# planepos = "<< planepos << " mm" << endl
    << "# panelangle = "<< panelangle << " rad" << endl
    << "# panelpos = "<< panelpos << " mm" << endl
    << "# wiredV = "<< wiredV << " mm" << endl
    << "# wiredW = "<< wiredW << " mm" << endl
    << "# twist = " << twist << " rad" << endl
    << "# xskew = "<< xskew << " rad" << endl
    << "# yskew = "<< yskew << " rad" << endl
    << "# zsqueeze = "<< zsqueeze << " mm" << endl
    << "# rsqueeze = "<< rsqueeze << " mm"  << endl
    << endl;
  for(unsigned ibuf=0;ibuf<3;ibuf++)
    os << "#" << endl;

  // global tracker
  os << "TABLE TrkAlignTracker " << endl;
  os << "#row, strawid, dx, dy, dz, rx, ry, rz" << endl;
  double dalphax = myrand.Gaus(0.0,trackerangle);
  double dalphay = myrand.Gaus(0.0,trackerangle);
  double dalphaz = myrand.Gaus(0.0,trackerangle);
  double dx = myrand.Gaus(0.0,trackerpos);
  double dy = myrand.Gaus(0.0,trackerpos);
  double dz = myrand.Gaus(0.0,trackerpos);
  os << "0, 0_0_0, "
    << dx << ", "
    << dy << ", "
    << dz << ", "
    << dalphax << ", "
    << dalphay << ", "
    << dalphaz << endl;

  for(unsigned ibuf=0;ibuf<3;ibuf++)
    os << "#" << endl;

  os << "TABLE TrkAlignPlane" << endl;
  os << "#row, strawid, dx, dy, dz, rx, ry, rz" << endl;
  for(unsigned iplane=0;iplane < NPlanes; ++iplane){
    double planez = (iplane-float(NPlanes)/2.0)*planezgap;
    double zfactor = planez/trackerlen; // scale by z position
    double dalphax = myrand.Gaus(xskew,planeangle);
    double dalphay = myrand.Gaus(yskew,planeangle);
    double dalphaz = myrand.Gaus(twist*zfactor,planeangle);
    double dx = myrand.Gaus(xskew*zfactor,planepos);
    double dy = myrand.Gaus(yskew*zfactor,planepos);
    double dz = myrand.Gaus(zsqueeze*zfactor,planepos);

    os << iplane << ", " << iplane << "_0_0, "
      << dx << ", "
      << dy << ", "
      << dz << ", "
      << dalphax << ", "
      << dalphay << ", "
      << dalphaz << endl;
  }
  for(unsigned ibuf=0;ibuf<3;ibuf++)
    os << "#" << endl;

  os << "TABLE TrkAlignPanel" << endl;
  os << "#row, strawid, dU, dV, dW, rU, rV, rW" << endl;
  int index(0);
  for(unsigned iplane=0;iplane < NPlanes; ++iplane){
    for(unsigned ipanel=0;ipanel < NPanels; ++ipanel){
      double dalphaU = myrand.Gaus(0,panelangle);
      double dalphaV = myrand.Gaus(0,panelangle);
      double dalphaW = myrand.Gaus(0,panelangle);
      double dU = myrand.Gaus(0.0,panelpos);
      double dV = myrand.Gaus(rsqueeze,panelpos);
      double dW = myrand.Gaus(0.0,panelpos);

      os << index << ", " << iplane << "_" << ipanel << "_0, "
        << dU << ", "
        << dV << ", "
        << dW << ", "
        << dalphax << ", "
        << dalphay << ", "
        << dalphaz << endl;
       index++;
   }
  }
  for(unsigned ibuf=0;ibuf<3;ibuf++)
    os << "#" << endl;
  os << "TABLE TrkAlignStraw" << endl;
  os << "#index,StrawId,wire_cal_dV,wire_cal_dW,wire_hv_dV,wire_hv_dW,straw_cal_dV,straw_cal_dW,straw_hv_dV,straw_hv_dW" << endl;
  index = 0;
  for(unsigned iplane=0;iplane < NPlanes; ++iplane){
    for(unsigned ipanel=0;ipanel<NPanels; ++ipanel){
      for(unsigned istraw=0;istraw<NStraws; ++istraw){
        double wire_cal_dV = myrand.Gaus(0.0,wiredV);
        double wire_cal_dW= myrand.Gaus(0.0,wiredW);
        double wire_hv_dV = myrand.Gaus(0.0,wiredV);
        double wire_hv_dW= myrand.Gaus(0.0,wiredW);
        double straw_cal_dV = myrand.Gaus(0.0,wiredV);
        double straw_cal_dW= myrand.Gaus(0.0,wiredW);
        double straw_hv_dV = myrand.Gaus(0.0,wiredV);
        double straw_hv_dW= myrand.Gaus(0.0,wiredW);
        os << index << ", " << iplane << "_" << ipanel << "_" << istraw << ", "
          << wire_cal_dV << ", " << wire_cal_dW << ", " << wire_hv_dV << ", " << wire_hv_dW << ", "
          << straw_cal_dV << ", " << straw_cal_dW << ", " << straw_hv_dV << ", " << straw_hv_dW << endl;
        index++;
      }
    }
  }
}
