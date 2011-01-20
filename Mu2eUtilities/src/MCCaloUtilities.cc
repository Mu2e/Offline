//
// $Id: MCCaloUtilities.cc,v 1.1 2011/01/20 21:28:39 onoratog Exp $
// $Author: onoratog $
// $Date: 2011/01/20 21:28:39 $
//
// Original author Gianni Onorato
//

// C++ includes
#include<iostream>

// Framework includes
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"

// Mu2e includes
#include "Mu2eUtilities/inc/MCCaloUtilities.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "ToyDP/inc/PhysicalVolumeInfoCollection.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "ToyDP/inc/StatusG4.hh"


using namespace std;

namespace mu2e {

  MCCaloUtilities::MCCaloUtilities()
  {
    _localRO = 0;
    _localCrystal = 0;
    _localVane = 0;
  }

  MCCaloUtilities::~MCCaloUtilities()
  {
  }


  void MCCaloUtilities::printOutCaloInfo() {

    edm::Service<GeometryService> geom;
    GeomHandle<Calorimeter> cg;
    
    double Hsize = cg->crystalHalfSize();
    double Hleng = cg->crystalHalfLength();
    
    double ROsize = cg->roHalfSize();
    
    cout << "Crystal HSize " << Hsize
         << "\nCrystal HLeng " << Hleng 
         << "\nRO size " << ROsize << endl;
    
    int nRO = cg->nRO();
    
    double minc[4][3], maxc[4][3];
    
    for (int i=0; i<4; ++i) {
      for (int k=0; k<4; ++k) {
        minc[i][k] = 100000;
        maxc[i][k] = -100000;
      }
    }
    
    for (int j=0; j< nRO/2; ++j) {
      
      int thevane = cg->getVaneByRO(2*j);
      
      cout << "Crystal n. " << cg->getCrystalByRO(2*j);
      cout << "\tVane " << thevane;
      
      
      CLHEP::Hep3Vector cntr = cg->getCrystalOriginByRO(2*j);
      CLHEP::Hep3Vector Xaxis = cg->getCrystalAxisByRO(2*j);
      CLHEP::Hep3Vector Yaxis = cg->getCrystalAxisByRO(2*j).orthogonal();
      CLHEP::Hep3Vector Zaxis(0,0,1);
      
      cout << "\tcenter " << cntr
           << "\tXaxis " << Xaxis 
           << "\tYaxis " << Yaxis;
           
      CLHEP::Hep3Vector toLeftC = ( (-Hsize) * Yaxis ) + ( (-Hleng) * Xaxis ) + ( (-Hsize) * Zaxis );

      CLHEP::Hep3Vector toRightC = ( (Hsize) * Yaxis ) + ( (Hleng+(2*ROsize)) * Xaxis ) + ( (Hsize) * Zaxis );

      cout << "toleft" << toLeftC
           << "\ttoright" << toRightC;

      
      CLHEP::Hep3Vector lcorner = cntr + toLeftC;
      CLHEP::Hep3Vector rcorner = cntr + toRightC;
      
      cout << "\tlcrn " << lcorner << "\trcrn " << rcorner << endl;
      
      if (minc[thevane][0] > lcorner.getX()) minc[thevane][0] = lcorner.getX();
      if (minc[thevane][0] > rcorner.getX()) minc[thevane][0] = rcorner.getX();
      if (minc[thevane][1] > lcorner.getY()) minc[thevane][1] = lcorner.getY();
      if (minc[thevane][1] > rcorner.getY()) minc[thevane][1] = rcorner.getY();
      if (minc[thevane][2] > lcorner.getZ()) minc[thevane][2] = lcorner.getZ();
      if (minc[thevane][2] > rcorner.getZ()) minc[thevane][2] = rcorner.getZ();
      
      if (maxc[thevane][0] < lcorner.getX()) maxc[thevane][0] = lcorner.getX();
      if (maxc[thevane][0] < rcorner.getX()) maxc[thevane][0] = rcorner.getX();
      if (maxc[thevane][1] < lcorner.getY()) maxc[thevane][1] = lcorner.getY();
      if (maxc[thevane][1] < rcorner.getY()) maxc[thevane][1] = rcorner.getY();
      if (maxc[thevane][2] < lcorner.getZ()) maxc[thevane][2] = lcorner.getZ();
      if (maxc[thevane][2] < rcorner.getZ()) maxc[thevane][2] = rcorner.getZ();
      
    } 
    
    cout << "Ranges: " << endl;
    for (int i=0; i<4; ++i) {
      cout << "\nVane " << i << " , min ( ";
      for (int k=0; k<3; ++k) {
        cout << minc[i][k] << (k==2?" ) ":", ");
      }
      cout << "   max ( ";
      for (int k=0; k<3; ++k) {
        cout << maxc[i][k] << (k==2?" ) ":", ");
      }
    }
    cout << endl;  
  }
  
  void MCCaloUtilities::setTrackAndRO(const edm::Event & event,
                                      SimParticleCollection::key_type track,
                                      uint32_t RO){

    _localRO = RO;
    edm::Service<GeometryService> geom;
    GeomHandle<Calorimeter> cg;
    _localCrystal  = cg->getCrystalByRO(_localRO);
    _localVane = cg->getVaneByRO(_localRO);
    
    edm::Handle<SimParticleCollection> simParticles;
    event.getByType(simParticles);

    SimParticle const& sim = simParticles->at(track);

    CLHEP::Hep3Vector origin = sim.startPosition();

    _startingVane = getStartingVane(origin);

    _fromOutside = (_startingVane == -1);
    
    _primary = !sim.hasParent();

    _generated = sim.fromGenerator();

    //    edm::Handle<PhysicalVolumeInfoCollection> volumes;
    //event.getRun().getByType(volumes);

    //    PhysicalVolumeInfo const& volInfob = volumes->at(sim.startVolumeIndex());
    //PhysicalVolumeInfo const& volInfoe = volumes->at(sim.endVolumeIndex());
    //    cout << "start: " << sim.startVolumeIndex() << "   " << volInfob.name() << "  " << volInfob.copyNo() << endl; 
    //cout << "end:   " << sim.endVolumeIndex() << "    " << volInfoe.name() << "  " << volInfoe.copyNo() << endl; 
    //  cout << "Start position: " << sim.startPosition() << '\n'
    //     << "End position:   " << sim.endPosition() << endl;
    //cout << "Particle process code " << sim.creationCode().name() << endl;
    //bool ID(false); 
    // if ( sim.hasParent()) {
    //   cout  << "Parent id " << sim.parentId() << endl;
    //   ID = true;
    // }  
    // if ( sim.fromGenerator() ) {
    //   cout << "Is from generator" << endl;
    //   ID =  true;
    // } 
    // if  (!ID) {
    //   cout << "dunno where it come from" << endl;
    // } 
     
  } 

  bool MCCaloUtilities::fromOutside() {
    return _fromOutside;
  }

  bool MCCaloUtilities::primary() {
    return _primary;
  }

  bool MCCaloUtilities::generated() {
    return _generated;
  }

  int MCCaloUtilities::startingVane() {
    return _startingVane;
  }

  int MCCaloUtilities::localVane() {
    return _localVane;  
  }

  int MCCaloUtilities::getStartingVane(CLHEP::Hep3Vector origin) {

    if (origin.getZ() >= 11950 && origin.getZ() <= 13450) {

      if (origin.getX() >= -4564 && origin.getX() <= -4264) {
        if (origin.getY() >= -55 && origin.getY() <= 55) {
          return 0;
        }
      }

      if (origin.getX() >= -3959 && origin.getX() <= -3849) {
        if (origin.getY() >= -660 && origin.getY() <= -360) {
          return 1;
        }
      }

      if (origin.getX() >= -3544 && origin.getX() <= -3244) {
        if (origin.getY() >= -55 && origin.getY() <= 55) {
          return 2;
        }
      }

      if (origin.getX() >= -3959 && origin.getX() <= -3849) {
        if (origin.getY() >= 360 && origin.getY() <= 660) {
          return 3;
        }
      }
    } 
    
    return -1;
  
  }

} // end namespace mu2e
