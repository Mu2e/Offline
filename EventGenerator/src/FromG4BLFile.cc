//
// Read particles from a file in G4beamline input format.
//
// $Id: FromG4BLFile.cc,v 1.7 2010/12/11 04:50:10 logash Exp $
// $Author: logash $ 
// $Date: 2010/12/11 04:50:10 $
//
// Original author Rob Kutschke
//
// The position is given in the Mu2e coordinate system.
// 

// Zlimits: 1685.84 1843.21
//   Length = 157.37           ( halflength=80)
//   Mid    = 1764.525

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
#include "GeometryService/inc/GeometryService.hh"
#include "Mu2eUtilities/inc/PDGCode.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "ToyDP/inc/G4BeamlineInfoCollection.hh"

// Root includes
#include "TH1F.h"
#include "TNtuple.h"

// Other external includes.
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"

using namespace std;

namespace mu2e {

  FromG4BLFile::FromG4BLFile( edm::Run const& , const SimpleConfig& config ):

    // Base class.
    GeneratorBase(),

    // From run time configuration file.
    _mean(config.getDouble("fromG4BLFile.mean",-1.)),
    _inputFileName(config.getString("fromG4BLFile.filename")),
    _doHistograms(config.getBool("fromG4BLFile.doHistograms", false)),
    _targetFrame(config.getBool("fromG4BLFile.targetFrame", true)),

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

    // Sanity check.
    if ( std::abs(_mean) > 99999. ) {
      throw cms::Exception("RANGE")
        << "FromG4BLFile has been asked to produce a crazily large number of particles.\n";
    }
    
    if ( config.hasName("fromG4BLFile.offset") ) {
      _offset = config.getHep3Vector("fromG4BLFile.offset");
    } else {
      _offset = CLHEP::Hep3Vector(0.0,0.0,1764.5);
    }

    // This should really come from the geometry service, not directly from the config file.
    // Or we should change this code so that its reference point is the production target midpoint
    edm::Service<GeometryService> geom;
    SimpleConfig const& geomConfig = geom->config();
    _prodTargetCenter = geomConfig.getHep3Vector("productionTarget.position");

    // Book histograms if enabled.
    if ( !_doHistograms ) return;

    edm::Service<edm::TFileService> tfs;

    edm::TFileDirectory tfdir = tfs->mkdir( "FromG4BLFile" );
    _hMultiplicity = tfdir.make<TH1F>( "hMultiplicity", "From G4BL file: Multiplicity",    20,  0.,  20.);

    _hMomentum     = tfdir.make<TH1F>( "hMomentum",     "From G4BL file: Momentum (MeV)",  100, 0., 1000. );

    _hCz           = tfdir.make<TH1F>( "hCz", "From G4BL file: cos(theta)",  100, -1.,  1.);

    _hX0           = tfdir.make<TH1F>( "hX0", "From G4BL file: X0",   100,  -40.,  40.);
    _hY0           = tfdir.make<TH1F>( "hY0", "From G4BL file: Y0",   100,  -20.,  20.);
    _hZ0           = tfdir.make<TH1F>( "hZ0", "From G4BL file: Z0",   100, -100., 100.);
    _hT0           = tfdir.make<TH1F>( "hT0", "From G4BL file: Time", 100, -500., 500.);

    _ntup          = tfdir.make<TNtuple>( "ntup", "G4BL Track ntuple",
                                          "x:y:z:p:cz:phi:pt:t:id:evtId:trkID:ParId:w");

  }

  FromG4BLFile::~FromG4BLFile(){
  }

  void FromG4BLFile::generate(ToyGenParticleCollection& genParts) {
    generate(genParts,0);
  }

  void FromG4BLFile::generate( ToyGenParticleCollection& genParts, G4BeamlineInfoCollection *extra ){

    // How many tracks in this event?
    long n = _mean < 0 ? static_cast<long>(-_mean): _randPoissonQ.fire();
    if ( _doHistograms ){
      _hMultiplicity->Fill(n);
    }

    // Particle data table.
    ConditionsHandle<ParticleDataTable> pdt("ignored");

    // Ntuple buffer.
    float nt[_ntup->GetNvar()];

    // Loop over all of the requested particles.
    for ( int j =0; j<n; ++j ){

      // Format of one line is: x y z Px Py Pz t PDGid EventID TrackID ParentID Weight
      double x, y, z, px, py, pz, t,  weight;
      int id, evtid, trkid, parentid;
      _inputFile >> x >> y >> z >> px >> py >> pz >> t >> id >> evtid >> trkid >> parentid >> weight;
      if ( !_inputFile ){
        throw cms::Exception("EOF")
          << "FromG4BLFile has reached an unexpected end of flie.\n";
      }

      // Express pdgId as the correct type.
      PDGCode::type pdgId = static_cast<PDGCode::type>(id);

      // 3D position in Mu2e coordinate system.
      CLHEP::Hep3Vector pos(x,y,z);
      if( _targetFrame ) {
	pos -= _offset;           // Move to target coordinate system
	pos += _prodTargetCenter; // Move to Mu2e coordinate system
      }

      // 4 Momentum.
      double mass = pdt->particle(id).mass().value();
      double e    = sqrt( px*px + py*py + pz*pz + mass*mass);
      CLHEP::HepLorentzVector p4(px,py,pz,e);

      std::cout <<x<<" "<<y<<" "<<z<<" "<<px<<" "<<py<<" "<<pz<<" "<<mass<<" "<<e<<endl;

      // Add particle to the output collection.
      genParts.push_back( ToyGenParticle( pdgId, GenId::fromG4BLFile, pos, p4, t) );

      // Add extra information to the output collection.
      if( extra ) {
	extra->push_back( G4BeamlineInfo(evtid,trkid,weight,t) );
      }

      if ( _doHistograms ) {

        // Magnitude of momentum, cos(theta) and azimuth.
        double p   = p4.vect().mag();
        double pt  = p4.vect().perp();
        double cz  = p4.vect().cosTheta();
        double phi = p4.vect().phi();

        _hMomentum->Fill(p);
        _hCz->Fill( cz );
        _hX0->Fill( x );
        _hY0->Fill( y );
        _hZ0->Fill( z );
        _hT0->Fill( t );

        nt[0]  = x;
        nt[1]  = y;
        nt[2]  = z;
        nt[3]  = p;
        nt[4]  = cz;
        nt[5]  = phi;
        nt[6]  = pt;
        nt[7]  = t;
        nt[8]  = id;
        nt[9]  = evtid;
        nt[10] = trkid;
        nt[11] = parentid;
        nt[12] = weight;
        _ntup->Fill(nt);
      }

    }

  } // end of generate

}
