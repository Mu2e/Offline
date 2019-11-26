#include "TRandom3.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>

using namespace std;
void RandomDeadStraws(unsigned ndeadplane, unsigned ndeadpanel, unsigned ndeadstraw, const char* outfile) {
  unsigned NStraws(96);
  unsigned NPlanes(36);
  unsigned NPanels(6);
  vector<string> deadstraws;
  vector<unsigned> inactiveplanes;
  vector<string> inactivepanels;
  vector<string> inactivestraws;
  
  TRandom3 myrand(3499803);

  filebuf fb;
  fb.open (outfile,ios::out);
  ostream os(&fb);
  // write header
  os << "physics.producers.makeSD.deadStrawList : {" << endl;

  if(ndeadplane > 0){
    os  << "deadPlane : [ " << endl;
    unsigned ndplane(0); 
    while(ndplane < ndeadplane){
      int plane = myrand.Integer(NPlanes);
      if(find(inactiveplanes.begin(), inactiveplanes.end(), plane) == inactiveplanes.end()){
	++ndplane;
	os  << "\"" << plane << "\"";
	if(ndplane < ndeadplane) os << "," << endl;
	inactiveplanes.push_back(plane);
      }
    }
    os  << "]" << endl;
  }

  if(ndeadpanel > 0){
    os  << "deadPanels : [ " << endl;
    unsigned ndpanel(0); 
    while(ndpanel < ndeadpanel){
      int plane = myrand.Integer(NPlanes);
      if(find(inactiveplanes.begin(), inactiveplanes.end(), plane) == inactiveplanes.end()){
	unsigned panel = myrand.Integer(NPanels);
	char panelid[80];
	snprintf(panelid,80,"\"%u_%u\"",plane,panel);
	if(find(inactivepanels.begin(),inactivepanels.end(),panelid) == inactivepanels.end()){
	  ++ndpanel;
	  os  << panelid;
	  if(ndpanel < ndeadpanel) os << "," << endl;
	  inactivepanels.push_back(panelid);
	}
      }
    }
    os  << "]" << endl;
  }

  if(ndeadstraw > 0) {
    os  << "deadStraws : [ " << endl;
    unsigned ndstraw(0);
    while (ndstraw < ndeadstraw){ 
      int plane = myrand.Integer(NPlanes);
      if(find(inactiveplanes.begin(), inactiveplanes.end(), plane) == inactiveplanes.end()){
     	unsigned panel = myrand.Integer(NPanels);
	char panelid[80];
	snprintf(panelid,80,"\"%u_%u\"",plane,panel);
	if(find(inactivepanels.begin(),inactivepanels.end(),panelid) == inactivepanels.end()){
	  unsigned straw = myrand.Integer(NStraws);
	  char strawid[80];
	  snprintf(strawid,80,"\"%u_%u_%u\"",plane,panel,straw);
	  if(find(deadstraws.begin(), deadstraws.end(), strawid) == deadstraws.end()){
	    ++ndstraw;
	    os << strawid;
	    if(ndstraw < ndeadstraw) os << "," << endl;
	    deadstraws.push_back(strawid);
	  }
	}
      }
    }
    os  << "]" << endl;;
  }
}

