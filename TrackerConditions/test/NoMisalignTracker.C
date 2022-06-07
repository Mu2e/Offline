#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
// generate parameters to define null alignment
// David Brown (LBNL)
//

using namespace std;
void NoMisalignTracker(const char* outfile="NoMisalignTracker.txt") {
  unsigned NPlanes(36);
  unsigned NPanels(6);
  unsigned NStraws(96);

  filebuf fb;
  fb.open (outfile,ios::out);
  ostream os(&fb);
  // write header
  os << "# Tracker With No Misalignments" << endl;
  for(unsigned ibuf=0;ibuf<3;ibuf++)
    os << "#" << endl;
  os << " TABLE TrkAlignTracker " << endl;
  os << "#row, strawid, dx, dy, dz, rx, ry, rz" << endl;
  os << "0, 0_0_0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0" << endl;
  for(unsigned ibuf=0;ibuf<3;ibuf++)
    os << "#" << endl;

  os << "TABLE TrkAlignPlane" << endl;
  os << "#row, strawid, dx, dy, dz, rx, ry, rz" << endl;
  for(unsigned iplane=0;iplane < NPlanes; ++iplane){
    os << iplane << ", " << iplane << "_0_0, "
      << "0.0, 0.0, 0.0, 0.0, 0.0, 0.0" << endl;
  }

  for(unsigned ibuf=0;ibuf<3;ibuf++)
    os << "#" << endl;
  os << "TABLE TrkAlignPanel" << endl;
  os << "#row, strawid, dU, dV, dW, rU, rV, rW" << endl;
  for(unsigned iplane=0;iplane < NPlanes; ++iplane){
    for(unsigned ipanel=0;ipanel<NPanels; ++ipanel){
      unsigned irow = iplane*NPanels + ipanel;
      os << irow << ", " << iplane << "_" << ipanel << "_0, "
        << "0.0, 0.0, 0.0, 0.0, 0.0, 0.0" << endl;
    }
  }

  for(unsigned ibuf=0;ibuf<3;ibuf++)
    os << "#" << endl;
  os << "TABLE TrkAlignStraw" << endl;
  os << "#index,StrawId,wire_cal_dV,wire_cal_dW,wire_hv_dV,wire_hv_dW,straw_cal_dV,straw_cal_dW,straw_hv_dV,straw_hv_dW" << endl;
  for(unsigned iplane=0;iplane < NPlanes; ++iplane){
    for(unsigned ipanel=0;ipanel<NPanels; ++ipanel){
      for(unsigned istraw=0;istraw<NStraws; ++istraw){
        unsigned irow = (iplane*NPanels + ipanel)*NStraws + istraw;
        os << irow << ", " << iplane << "_" << ipanel << "_" << istraw << ", "
          << "0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0" << endl;
      }
    }
  }
}
