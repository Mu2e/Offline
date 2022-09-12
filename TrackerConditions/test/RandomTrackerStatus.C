//
//  Generate a randomized TrackerStatus text file given a number of dead straws, plans, and panels, a
//  give number of noisy channels.  Invoke this as (for example) :
//  root> .L Offline/TrackerConditions/test/RandomTrackerStatus.C
//  root> RandomTrackerStatus(102942,1,3,288,288,"RandomTrackerStatus.txt")
//
//  Then use it in an Offline reconstruction (or simulation digitization) job as
//
#include "TRandom3.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>

using namespace std;
void RandomTrackerStatus(unsigned seed, unsigned ndeadplanes, unsigned ndeadpanels, unsigned ndeadstraws, unsigned nnoisystraws, const char* outfile="RandomTrackerStatus.txt") {
  unsigned NStraws(96); // straws/plane
  unsigned NPanels(6); // panels/plane
  unsigned NPlanes(36);  // planes
  vector<unsigned> deadplanes;
  vector<unsigned> deadpanels;
  vector<unsigned> deadstraws;
  vector<unsigned> noisystraws;

  TRandom3 myrand(seed);

  filebuf fb;
  fb.open (outfile,ios::out);
  ostream os(&fb);
  // write header
  os << "#" << endl;
  os << "# random TrackerStatus generated with " << ndeadplanes << " dead planes, "
    << ndeadpanels << " dead panels, "
    << ndeadstraws << " dead straws, "
    << nnoisystraws << " noisy straws" << endl;
  os << "#" << endl;

  os  << "TABLE TrkPlaneStatus" << endl;
  if(ndeadplanes > 0){
    unsigned ndplane(0);
    while(ndplane < ndeadplanes){
      int plane = myrand.Integer(NPlanes);
      if(find(deadplanes.begin(), deadplanes.end(), plane) == deadplanes.end()){
        ++ndplane;
        os  << std::setw(2) << plane << "_0_0, Absent" << std::endl;
        deadplanes.push_back(plane);
      }
    }
  }
  os << "#" << endl << "#" << endl << "#" << endl;

  os  << "TABLE TrkPanelStatus"  << endl;
  if(ndeadpanels > 0){
    unsigned ndpanel(0);
    while(ndpanel < ndeadpanels){
      int plane = myrand.Integer(NPlanes);
      if(find(deadplanes.begin(), deadplanes.end(), plane) == deadplanes.end()){
        unsigned panel = myrand.Integer(NPanels);
        unsigned upanel = plane*NPanels + panel;
        if(find(deadpanels.begin(),deadpanels.end(),upanel) == deadpanels.end()){
          ++ndpanel;
          char panelid[80];
          snprintf(panelid,80,"%u_%u_0, NoHV",plane,panel);
          os  << panelid << std::endl;
          deadpanels.push_back(upanel);
        }
      }
    }
  }
  os << "#" << endl << "#" << endl << "#" << endl;

  os  << "TABLE TrkStrawStatusLong" << endl;
  if(ndeadstraws > 0) {
    unsigned ndstraw(0);
    while (ndstraw < ndeadstraws){
      int plane = myrand.Integer(NPlanes);
      if(find(deadplanes.begin(), deadplanes.end(), plane) == deadplanes.end()){
        unsigned panel = myrand.Integer(NPanels);
        unsigned upanel = plane*NPanels + panel;
        if(find(deadpanels.begin(),deadpanels.end(),upanel) == deadpanels.end()){
          unsigned straw = myrand.Integer(NStraws);
          unsigned ustraw = plane*NStraws*NPanels + panel*NStraws + straw;
          if(find(deadstraws.begin(), deadstraws.end(), ustraw) == deadstraws.end()){
            ++ndstraw;
            char strawid[80];
            snprintf(strawid,80,"%u_%u_%u, NoWire",plane,panel,straw);
            os << strawid << std::endl;
            deadstraws.push_back(ustraw);
          }
        }
      }
    }
  }
  os << "#" << endl << "#" << endl << "#" << endl;

  os  << "TABLE TrkStrawStatusShort" << endl;
  if(nnoisystraws > 0) {
    unsigned ndstraw(0);
    while (ndstraw < nnoisystraws){
      int plane = myrand.Integer(NPlanes);
      if(find(deadplanes.begin(), deadplanes.end(), plane) == deadplanes.end()){
        unsigned panel = myrand.Integer(NPanels);
        unsigned upanel = plane*NPanels + panel;
        if(find(deadpanels.begin(),deadpanels.end(),upanel) == deadpanels.end()){
          unsigned straw = myrand.Integer(NStraws);
          unsigned ustraw = plane*NStraws*NPanels + panel*NStraws + straw;
          if(find(noisystraws.begin(), noisystraws.end(), ustraw) == noisystraws.end()
              && find(deadstraws.begin(), deadstraws.end(), ustraw) == deadstraws.end()){
            ++ndstraw;
            char strawid[80];
            snprintf(strawid,80,"%u_%u_%u, Noise",plane,panel,straw);
            os << strawid << std::endl;
            noisystraws.push_back(ustraw);
          }
        }
      }
    }
  }
}
