//
// Read particles from a file in G4beamline input format.
//
// $Id: FromG4BLFile.cc,v 1.1 2010/08/30 22:50:00 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/08/30 22:50:00 $
//
// Original author Rob Kutschke
//
// The position is given in the Mu2e coordinate system.
// 

#include <iostream>

// Framework includes
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Services/interface/TFileService.h"
#include "FWCore/Framework/interface/TFileDirectory.h"

// Mu2e includes
#include "EventGenerator/inc/FromG4BLFile.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "LTrackerGeom/inc/LTracker.hh"
#include "ITrackerGeom/inc/ITracker.hh"
#include "Mu2eUtilities/inc/PDGCode.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"

// Root includes
#include "TH1F.h"

// Other external includes.
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"

using namespace std;

//using CLHEP::Hep3Vector;
//using CLHEP::HepLorentzVector;
//using CLHEP::RandFlat;
//using CLHEP::twopi;


namespace mu2e {

  FromG4BLFile::FromG4BLFile( edm::Run const& , const SimpleConfig& config ):

    // Base class.
    GeneratorBase(),

    // From run time configuration file.
    _mean(config.getDouble("particleGun.mean",-1.)),
    _inputFileName(config.getString("fromG4BLFile.filename")),
    _doHistograms(config.getBool("particleGun.doHistograms", false)),

    // Random number distributions; getEngine() comes from base class.
    _randPoissonQ( getEngine(), std::abs(_mean) ),

    // Open the input file.
    _inputFile(_inputFileName.c_str()),

    // Histogram pointers
    _hMultiplicity(0),
    _hMomentum(0),
    _hCz(0),
    _hX0(0),
    _hY0(0),
    _hZ0(0),
    _hT0(0){
    
    // Book histograms if enabled.
    if ( !_doHistograms ) return;

    edm::Service<edm::TFileService> tfs;

    edm::TFileDirectory tfdir = tfs->mkdir( "FromG4BLFile" );
    _hMultiplicity = tfdir.make<TH1F>( "hMultiplicity", "From G4BL file: Multiplicity",    20,  0.,  20.);

    _hMomentum     = tfdir.make<TH1F>( "hMomentum",     "From G4BL file: Momentum (MeV)",  100, 0., 1000. );

    _hCz           = tfdir.make<TH1F>( "hCz",           "From G4BL file: cos(theta)",      100, -1.,  1.);

    _hX0           = tfdir.make<TH1F>( "hX0", "From G4BL file: X0",              100,  -4000.,    4000.);
    _hY0           = tfdir.make<TH1F>( "hY0", "From G4BL file: Y0",              100,  -1000.,    1000.);
    _hZ0           = tfdir.make<TH1F>( "hZ0", "From G4BL file: Z0",              100, -10000.,   10000.);
    _hT0           = tfdir.make<TH1F>( "hT0", "From G4BL file: Time",            100,      0.,    1800.);

  }

  FromG4BLFile::~FromG4BLFile(){
  }

  void FromG4BLFile::generate( ToyGenParticleCollection& genParts ){

    long n = _mean < 0 ? static_cast<long>(-_mean): _randPoissonQ.fire();
    if ( _doHistograms ){
      _hMultiplicity->Fill(n);
    }

    for ( int j =0; j<n; ++j ){

      CLHEP::Hep3Vector pos;

      // 4 Momentum.
      CLHEP::HepLorentzVector p4;

      // Time
      //double time = _tmin + _dt*_randFlat.fire();

      //genParts.push_back( ToyGenParticle( _pdgId, GenId::particleGun, pos, p4, time));

      /*
      cout << "Generated position: " 
           << pos << " "
           << p4 << " "
           << p4.vect().mag() << " "
           << time
           << endl;

      if ( _doHistograms ) {
        _hMomentum->Fill(p);
        _hCz->Fill( p4.vect().cosTheta());
        _hX0->Fill( pos.x() );
        _hY0->Fill( pos.y() );
        _hZ0->Fill( pos.z() );
        _hT0->Fill( time );
      }
      */

    }

  } // end of generate

}
