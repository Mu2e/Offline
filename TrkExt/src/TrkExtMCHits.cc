//
//
//  $Id: TrkExtMCHits.cc,v 1.3 2013/05/16 18:23:39 mjlee Exp $
//  $Author: mjlee $
//  $Date: 2013/05/16 18:23:39 $
//
//  Original author MyeongJae Lee
//
//

// C++ includes.
#include <iostream>
#include <string>

// Framework includes.

#include "CLHEP/Vector/ThreeVector.h"
#include "TrkExt/inc/TrkExtMCHits.hh"
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "MCDataProducts/inc/SimParticle.hh"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"


using namespace CLHEP;
using namespace std;

namespace mu2e {

 
  TrkExtMCHits::TrkExtMCHits( ) { 
    _g4ModuleLabel = "";
    _instanceName = "";
    _hitcol.clear();
    _distcut = 1.;
  }

  TrkExtMCHits::~TrkExtMCHits() { }

  TrkExtMCHits::TrkExtMCHits(const art::Event& event, std::string g4ModuleLabel, std::string instanceName, cet::map_vector_key simid, double distcut) {
    _g4ModuleLabel = g4ModuleLabel;
    _instanceName = instanceName;
    _hitcol.clear();
    _distcut = distcut;

    art::Handle<StepPointMCCollection> hits;
    event.getByLabel(_g4ModuleLabel, _instanceName, hits);
    if (!hits.isValid()) {
      _hitcol.clear();
      cerr << "TrkExtMCHits : cannot find StepPointMC for " << instanceName << endl;
      return;
    }

    if (hits->size() <=0) {
      _hitcol.clear();
      return;
    }

    unsigned int i;
    for ( i=0; i<hits->size(); ++i ){
      const StepPointMC& hit = (*hits)[i];
      if (hit.trackId() == simid) {
        bool flag = false;
        for (unsigned int k = 0 ; k < _hitcol.size() ; ++k) {
          StepPointMCCollection & mccol = _hitcol[k];
          for (unsigned int l = 0 ; l < mccol.size() ; ++l) {
            StepPointMC & mc = mccol[l];
            if ((mc.position() - hit.position()).mag() <_distcut) {
              mccol.push_back(hit);
              flag = true;
              break;
            }
          }
          if (flag) break;
        }
        if (!flag) {
          StepPointMCCollection mccol;
          mccol.push_back(hit);
          _hitcol.push_back(mccol);
        }
      }
    }
    return;
  }

  Hep3Vector TrkExtMCHits::momentum (unsigned int clust) const {
    const vector<StepPointMC> & cluster = _hitcol[clust];
    if (cluster.size() <=0) {
      cerr << "TrkExtMCHits Error : invalid _hitcol" << endl;
      Hep3Vector p(0,0,0);
      return p;
    }
    else if (cluster.size() ==1) {
      return cluster[0].momentum();
    }
    else if (cluster.size() == 2) {
      double z1 = cluster[0].position().z();
      double z2 = cluster[1].position().z();
      double zc = 0.5 * (z1+z2);
      double px = interpolate2(zc, z1, z2, cluster[0].momentum().x(), cluster[1].momentum().x());
      double py = interpolate2(zc, z1, z2, cluster[0].momentum().y(), cluster[1].momentum().y());
      double pz = interpolate2(zc, z1, z2, cluster[0].momentum().z(), cluster[1].momentum().z());
      Hep3Vector p(px,py,pz);
      return p;
    }
    else {
      unsigned int i1 = 0;
      unsigned int i3 = cluster.size()-1;
      unsigned int i2 = (i1+i3)/2;
      double z1 = cluster[i1].position().z();
      double z2 = cluster[i2].position().z();
      double z3 = cluster[i3].position().z();
      double zc = (z1+z3)*0.5;
      if (z2 == zc) {
        return cluster[i2].momentum();
      }
      double px = interpolate3(zc, z1,z2,z3,cluster[i1].momentum().x(), cluster[i2].momentum().x(), cluster[i3].momentum().x());
      double py = interpolate3(zc, z1,z2,z3,cluster[i1].momentum().y(), cluster[i2].momentum().y(), cluster[i3].momentum().y());
      double pz = interpolate3(zc, z1,z2,z3,cluster[i1].momentum().z(), cluster[i2].momentum().z(), cluster[i3].momentum().z());
      Hep3Vector p(px,py,pz);
      return p;
    }
  }

  Hep3Vector TrkExtMCHits::position (unsigned int clust) const {
    const vector<StepPointMC> & cluster = _hitcol[clust];
    if (cluster.size() <=0) {
      cerr << "TrkExtMCHits Error : invalid _hitcol" << endl;
      Hep3Vector p(0,0,0);
      return p;
    }
    else if (cluster.size() ==1) {
      return cluster[0].position();
    }
    else if (cluster.size() == 2) {
      double z1 = cluster[0].position().z();
      double z2 = cluster[1].position().z();
      double zc = 0.5 * (z1+z2);
      double x = interpolate2(zc, z1, z2, cluster[0].position().x(), cluster[1].position().x());
      double y = interpolate2(zc, z1, z2, cluster[0].position().y(), cluster[1].position().y());
      Hep3Vector p(x,y,zc);
      return p;
    }
    else {
      unsigned int i1 = 0;
      unsigned int i3 = cluster.size()-1;
      unsigned int i2 = (i1+i3)/2;
      double z1 = cluster[i1].position().z();
      double z2 = cluster[i2].position().z();
      double z3 = cluster[i3].position().z();
      double zc = (z1+z3)*0.5;
      if (z2 == zc) {
        return cluster[i2].position();
      }
      double x = interpolate3(zc, z1,z2,z3,cluster[i1].position().x(), cluster[i2].position().x(), cluster[i3].position().x());
      double y = interpolate3(zc, z1,z2,z3,cluster[i1].position().y(), cluster[i2].position().y(), cluster[i3].position().y());
      Hep3Vector p(x,y,zc);
      return p;
    }
  }

  double TrkExtMCHits::time (unsigned int clust) const {
    const vector<StepPointMC> & cluster = _hitcol[clust];
    if (cluster.size() <=0) {
      cerr << "TrkExtMCHits Error : invalid _hitcol" << endl;
      return 0;
    }
    else {
      double ret = 0;
      for (unsigned int i = 0 ; i <cluster.size() ; ++i) {
        ret += cluster[i].time();
      }
      return ret / (double)(cluster.size());
    }
  }

  double TrkExtMCHits::deltap (unsigned int clust) const {
    const vector<StepPointMC> & cluster = _hitcol[clust];
    return fabs(cluster[0].momentum().mag() - cluster[cluster.size()-1].momentum().mag());
  }

  double TrkExtMCHits::eDep (unsigned int clust) const {
    const vector<StepPointMC> & cluster = _hitcol[clust];
    double eDep = 0;
    for (unsigned int i = 0 ; i < cluster.size() ; ++i) {
      eDep += cluster[i].eDep();
    }
    return eDep;
  }

  double TrkExtMCHits::ionizingEdep (unsigned int clust) const {
    const vector<StepPointMC> & cluster = _hitcol[clust];
    double eDep = 0;
    for (unsigned int i = 0 ; i < cluster.size() ; ++i) {
      eDep += cluster[i].ionizingEdep();
    }
    return eDep;
  }

  double TrkExtMCHits::nonIonizingEdep (unsigned int clust) const {
    const vector<StepPointMC> & cluster = _hitcol[clust];
    double eDep = 0;
    for (unsigned int i = 0 ; i < cluster.size() ; ++i) {
      eDep += cluster[i].nonIonizingEDep();
    }
    return eDep;
  }

  double TrkExtMCHits::interpolate3 (double z, double x1, double x2, double x3, double y1, double y2, double y3) const {
    double x1_x2 = x1 - x2;
    double x2_x3 = x2 - x3;
    double x3_x1 = x3 - x1;
    if (x1_x2 !=0 && x2_x3 != 0 && x3_x1 != 0) {
      return y1*(z-x2)*(x3-z)/x1_x2/x3_x1 + y2*(x1-z)*(z-x3)/x1_x2/x2_x3 + y3*(z-x1)*(x2-z)/x3_x1/x2_x3;
    }
    else {
      return 0;
    }
  }

  double TrkExtMCHits::interpolate2 (double z, double x1, double x2, double y1, double y2) const {
    if (x2 != x1) {
      return (y2-y1)/(x2-x1)*(z-x1)+y1;
    }
    else {
      return 0;
    }
  }

    
} // end namespace mu2e

