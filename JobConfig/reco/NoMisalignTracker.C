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
  os << "0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0" << endl;
  for(unsigned ibuf=0;ibuf<3;ibuf++)
    os << "#" << endl;
 
  os << "TABLE TrkAlignPlane" << endl;
  for(unsigned iplane=0;iplane < NPlanes; ++iplane){
     os << iplane << ", " 
       << "0.0, 0.0, 0.0, 0.0, 0.0, 0.0" << endl;
  }

  for(unsigned ibuf=0;ibuf<3;ibuf++)
    os << "#" << endl;
  os << "TABLE TrkAlignPanel" << endl;
  for(unsigned ipanel=0;ipanel<NPanels*NPlanes; ++ipanel){
    os << ipanel << ", "
      << "0.0, 0.0, 0.0, 0.0, 0.0, 0.0" << endl;
  }

  for(unsigned ibuf=0;ibuf<3;ibuf++)
    os << "#" << endl;
  os << "TABLE TrkAlignStraw" << endl;
  for(unsigned istraw=0;istraw<NStraws*NPanels*NPlanes; ++istraw){
    os << istraw << ", "
      << "0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0" << endl;
  }
}
