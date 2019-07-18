//
// Read particles from a file in G4beamline input format.
// Position of the GenParticles is in the Mu2e coordinate system.
//
// $Id: FromG4BLFile.cc,v 1.30 2014/01/15 17:07:30 tassiell Exp $
// $Author: tassiell $
// $Date: 2014/01/15 17:07:30 $
//
// Original author Rob Kutschke
//
// The position is given in the Mu2e coordinate system.
//
// Notes:
// 1) There are two modes that are selected by the config variable:
//      fromG4BLFile.targetFrame = true or false
//
//    True:
//       Particles should be created in the production target.
//       For historical reasons, some versions of the G4beamline
//       have the production target in different locations.
//       So the variable
//          fromG4BLFile.prodTargetOff
//       should be set to the position of the center of the production
//       target in the G4beamline coordinate system.  This is subtracted
//       from the position read from the file to yield a position that
//       is centered on the production target.  The code then uses
//       the GeometryService value of the target position to place the
//       generated particles in the event.
//
//    False:
//       Positions in the input file are specified in the G4beamline
//       coordinate system.  This is transformed to the Mu2e coordinate
//       system using the configuration parameter:
//              fromG4BLFile.g4beamlineOrigin
//       There is a second offset that can be used:
//              fromG4BLFile.g4beamlineExtraOffset
//       Normally this should be (0,0,0) but it can be used to fix
//       small deviations in geometry specifications between G4beamline
//       and Mu2esim.
//

#include <iostream>

// Framework includes
#include "art/Framework/Principal/Run.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "EventGenerator/inc/FromG4BLFile.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "MCDataProducts/inc/G4BeamlineInfoCollection.hh"
#include "ProductionTargetGeom/inc/ProductionTarget.hh"

// Root includes
#include "TH1F.h"
#include "TNtuple.h"

// Other external includes.
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

using namespace std;

namespace mu2e {

  FromG4BLFile::FromG4BLFile(CLHEP::HepRandomEngine& engine, art::Run const&, SimpleConfig const& config):

    // Base class.
    GeneratorBase(true),

    // From run time configuration file.
    _mean(config.getDouble("fromG4BLFile.mean",-1.)),
    _prodTargetOffset(),
    _prodTargetCenter(),
    _g4beamlineOrigin(),
    _g4beamlineExtraOffset(CLHEP::Hep3Vector()),
    _inputFileName(config.getString("fromG4BLFile.filename")),
    _pdgIdToKeep(),
    _doHistograms(config.getBool("fromG4BLFile.doHistograms", false)),
    _targetFrame(config.getBool("fromG4BLFile.targetFrame", false)),
    _throwOnOutOfWorld(config.getBool("fromG4BLFile.throwOnOutOfWorld", true)),
    _nPartToSkip(config.getInt("fromG4BLFile.particlesToSkip",0)),
    _duplicate(config.getBool("fromG4BLFile.duplicateParticles",false)),

    _randPoissonQ(engine, std::abs(_mean)),

    // Open the input file.
    //_inputFile(_inputFileName.c_str()),
    _inputFile(),

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
      throw cet::exception("RANGE")
        << "FromG4BLFile has been asked to produce a crazily large number of particles.\n";
    }

    CLHEP::Hep3Vector offset_default(0.0,0.0,1764.5);
    _prodTargetOffset = config.getHep3Vector("fromG4BLFile.prodTargetOff", offset_default);

    CLHEP::Hep3Vector g4beamlineOrigin_default(3904.,0.,-7929.);
    _g4beamlineOrigin = config.getHep3Vector("fromG4BLFile.g4beamlineOrigin", g4beamlineOrigin_default);

    CLHEP::Hep3Vector g4beamlineExtraOffset_default;
    _g4beamlineExtraOffset = config.getHep3Vector("fromG4BLFile.g4beamlineExtraOffset", g4beamlineExtraOffset_default);

    vector<int> default_pdgIdToKeep;
    config.getVectorInt("fromG4BLFile.pdgIdToKeep", _pdgIdToKeep, default_pdgIdToKeep );

    _prodTargetCenter = GeomHandle<ProductionTarget>()->position();

    // Construct the full path to the input file.  Accept absolute paths and paths
    // starting with "." as is; other wise apply the standard search path.
    string path(_inputFileName);
    if ( path.substr(0,1) != "/" && path.substr(0,1) != "." ){
      ConfigFileLookupPolicy configFile;
      path = configFile(_inputFileName);
    }

    // Open the input file.
    _inputFile.open(path.c_str());

    //Skip the first nParticlesToSkip of the file;
    //10000 is a fake number of chars, before reaching the newline character.
    for (int jid=0; jid < _nPartToSkip; ++jid) {
      _inputFile.ignore(10000, '\n');
    }

    // Book histograms if enabled.
    if ( !_doHistograms ) return;

    art::ServiceHandle<art::TFileService> tfs;

    art::TFileDirectory tfdir = tfs->mkdir( "FromG4BLFile" );
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

  void FromG4BLFile::generate(GenParticleCollection& genParts) {
    generate(genParts,0);
  }

  void FromG4BLFile::generate( GenParticleCollection& genParts, G4BeamlineInfoCollection *extra ){

    // How many tracks in this event?
    long n = _mean < 0 ? static_cast<long>(-_mean): _randPoissonQ.fire();
    if ( _doHistograms ){
      _hMultiplicity->Fill(n);
    }

    // Particle data table.
    GlobalConstantsHandle<ParticleDataTable> pdt;

    GeomHandle<WorldG4> worldG4;

    // Read particles from the file until the requested number of particle have been read.
    int j(0);
    while ( j < n ){

      if( _duplicate && j>0 ) break;

      // Format of one line from the input file is: x y z Px Py Pz t PDGid EventID TrackID ParentID Weight
      double x, y, z, px, py, pz, t,  weight;
      int id, evtid, trkid, parentid;

      // Invariant: pdgId of particle just read is not on the list of required pdgIds.
      //            and location of the particle is outside of the world volume.
      bool looking(true);
      while (looking){

        _inputFile >> x >> y >> z >> px >> py >> pz >> t >> id >> evtid >> trkid >> parentid >> weight;
        if ( !_inputFile ){
          throw cet::exception("EOF")
            << "FromG4BLFile has reached an unexpected end of file.\n";
        }

        // Check if this pdgId is on the allowed list; an empty list means keep everything.
        if ( !_pdgIdToKeep.empty() ) {

          bool idOK(false);
          for ( std::vector<int>::const_iterator i=_pdgIdToKeep.begin();
                i != _pdgIdToKeep.end(); ++i ){

            // Accept the particle just read.
            if ( id == *i ) {
              idOK = true;
              break;
            }
          }

          // This particle failed this test, so try the next particle.
          if ( !idOK ) continue;

        }

        // Check that the space point is inside the G4 world.
        bool inWorld(false);
        if ( worldG4->inWorld( CLHEP::Hep3Vector(x,y,z) ) ) {
          inWorld = true;
        } else{
          cout << "Not in the world, skipping: " << _throwOnOutOfWorld << endl;
          if ( _throwOnOutOfWorld ){
            throw cet::exception("GEOM")
              << "FromG4BLFile has read in a point that is outside of the G4 world: "
              << CLHEP::Hep3Vector(x,y,z)
              << "\n";
          }else{
            mf::LogInfo("GEOM")
              << "FromG4BLFile:  Skipping particle that is out of the world volume.  Position is "
              << CLHEP::Hep3Vector(x,y,z)
              << "\n";
          }
        }

        // This particle passes all tests so keep it; if not, try the next particle.
        if ( inWorld ) looking=false;

      } // end loop over input particles.
      ++j;

      // Express pdgId as the correct type.
      PDGCode::type pdgId = static_cast<PDGCode::type>(id);

      // 3D position in Mu2e coordinate system.
      CLHEP::Hep3Vector pos(x,y,z);
      CLHEP::Hep3Vector oldpos(pos);
      if( _targetFrame ) {
        pos -= _prodTargetOffset;           // Move to target coordinate system
        pos += _prodTargetCenter; // Move to Mu2e coordinate system
      } else{
        pos += _g4beamlineOrigin;
        pos += _g4beamlineExtraOffset;
      }

      // 4 Momentum.
      double mass = pdt->particle(id).ref().mass().value();
      double e    = sqrt( px*px + py*py + pz*pz + mass*mass);
      CLHEP::HepLorentzVector p4(px,py,pz,e);

      // Add particle to the output collection.
      for( int k=0; k<n; ++k ) {
        genParts.push_back( GenParticle( pdgId, GenId::fromG4BLFile, pos, p4, t) );
        if( !_duplicate ) break;
      }

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

        // Ntuple buffer.
        float nt[_ntup->GetNvar()];

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
