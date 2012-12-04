#ifndef CREATEWIRES_H
#define CREATEWIRES_H 1
#include "StrawHit.hh"
#include "StrawHitMCTruth.hh"
#include "IlcDCHcluster.h"
#include "TObjArray.h"
#include <vector>
using namespace std; 

void createWires();
void readData();
TObjArray *mu2eHits2ilc(const vector<mu2e::StrawHit> *shits,const vector<mu2e::StrawHitMCTruth> *shitmcs,const vector<int>* select=0,double t0=-10, bool useZCoordinate=false);

#endif
