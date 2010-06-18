//
// Module to understand how to use the BaBar Kalman filter package.
// Not for general use.
//
// $Id: KalmanT01_plugin.cc,v 1.4 2010/06/18 19:24:05 genser Exp $
// $Author: genser $
// $Date: 2010/06/18 19:24:05 $
//
// Original author Rob Kutschke
//
// Notes:
// 1) First version is designed to work with measurements only; ie
//    no multiple scattering.  It should work with the
//    

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>

// Framework includes.
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Services/interface/TFileService.h"
#include "FWCore/Framework/interface/TFileDirectory.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Root includes.
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TNtuple.h"

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "LTrackerGeom/inc/LTracker.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "ToyDP/inc/ToyGenParticleCollection.hh"
#include "ToyDP/inc/SimParticleCollection.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "Mu2eUtilities/inc/toHepPoint.hh"
#include "KalmanTests/inc/printTrkParams.hh"

// Babar Kalman filter includes
#include "TrkBase/TrkHelixUtils.hh"
#include "BField/BFieldFixed.hh"

// Other includes.
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Random/RandGaussQ.h"

using namespace std;
using edm::Event;
using CLHEP::Hep3Vector;


namespace mu2e {

  // Indices of track parametes within the 
  // For the () operator;
  // The [] operator uses 0 based numbers.  See TrkBase/HelixTraj.hh
  // Is there an official Babar version of this for 1 based indices?
  enum TrkParIdx { kd0 =1, kphi0, kom, kz0, kct };

  class KalmanT01 : public edm::EDAnalyzer {
  public:
    explicit KalmanT01(edm::ParameterSet const& pset) : 
      _generatorModuleLabel(pset.getUntrackedParameter<string>("generatorModuleLabel","g4run")),
      _diagLevel(pset.getUntrackedParameter<int>("diagLevel",1)),
      _maxFullPrint(pset.getUntrackedParameter<int>("maxFullPrint",5)),
      _hDriftDist(0),
      _hCheckPointRadius(0),
      _messageCategory("FitTest"){
    }
    virtual ~KalmanT01() { }

    virtual void beginJob(edm::EventSetup const&);
 
    void analyze( const edm::Event& e, edm::EventSetup const&);

  private:
    
    // Name of the module that made these hits.
    string _generatorModuleLabel;

    // Diagnostics level.
    int _diagLevel;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    // Pointers to diagnostic histograms.
    TH1F* _hDriftDist;
    TH1F* _hCheckPointRadius;

    TH1F* _hd0;
    TH1F* _hphi0;
    TH1F* _hk;
    TH1F* _hz0;
    TH1F* _hct;
    TH1F* _hs0;


    // A category for the error logger.
    const std::string _messageCategory;

  };

  void KalmanT01::beginJob(edm::EventSetup const& ){

    // Create histograms if diagnostics are enabled.
    if ( _diagLevel > 0 ){

      edm::Service<edm::TFileService> tfs;

      _hDriftDist = tfs->make<TH1F>( "hDriftDist", "Generated Drift Distance;(mm)", 100, -3., 3. );
      _hCheckPointRadius = tfs->make<TH1F>( "hCheckPointRadius",  "Radius of Reference point; (mm)",
                                            100, 2.4, 2.6 );
      _hd0   = tfs->make<TH1F>( "hd0",   "Gen d0;(cm)",    100, -10., 10. );

      _hphi0 = tfs->make<TH1F>( "hphi0", "Gen phi00;(radians)", 100, -CLHEP::pi, CLHEP::pi);
      _hz0   = tfs->make<TH1F>( "hz0",   "Gen z0;(cm)",     70, -500., -300. );
      _hct   = tfs->make<TH1F>( "hct",   "Gen ct;",     70, 0.3, 0.8 );
      _hs0   = tfs->make<TH1F>( "hs0",   "Gen s0;",     70, -10., 10. );
    }
  }

  void KalmanT01::analyze(const edm::Event& event, edm::EventSetup const&) {

    // Counter used by debug printout. 
    static int ncalls(0);
    ++ncalls;

    // Tracker origin in the Mu2e system.
    static CLHEP::Hep3Vector trackerOrigin( -3904., 0, 10200.);

    // Get a reference to one of the L or T trackers.  Throw exception if not successful.
    const Tracker& tracker = getTrackerOrThrow();

    // Ask the event to give us a handle to the requested hits.
    edm::Handle<StepPointMCCollection> points;
    event.getByLabel(_generatorModuleLabel,points);

    // Get handles to the generated and simulated particles.
    edm::Handle<ToyGenParticleCollection> gensHandle;
    event.getByType(gensHandle);
    const ToyGenParticleCollection& gens = *gensHandle;

    edm::Handle<SimParticleCollection> simsHandle;
    event.getByType(simsHandle);
    const SimParticleCollection& sims = *simsHandle;

    // Eventually get this from a calibration DB.
    static const double sigma=0.1*CLHEP::mm;

    const SimParticle& sim = sims.front();

    // COnvert to BaBar style helix params.
    // To test that I know how to use the tools.
    HepVector startPar(5);
    double s0;
    double charge(-1.);
    const BFieldFixed bfield(0.,0.,-1.,0.);
    const CLHEP::Hep3Vector pos0 = sim.startPosition() - trackerOrigin;

    TrkHelixUtils::helixFromMom( startPar, s0, 
                                 toHepPoint(pos0),
                                 sim.startMomentum(),
                                 charge,
                                 bfield);
    cout << "Starting position: "
         << pos0 << " "
         << endl;
    cout << "Track parms: ";
    printTrkParams(cout,startPar,false);
    cout << "Flight length: " << s0 << endl;
    _hd0->  Fill(startPar(kd0));
    _hphi0->Fill(startPar(kphi0));
    _hz0->Fill(startPar(kz0));
    _hct->Fill(startPar(kct));
    _hs0->Fill(s0);

    double ct = startPar(5);
    double sz = 1./sqrt( 1. + ct*ct );
    double cz = ct*cz;

    double radius = 1./std::abs(startPar(3));
    double arc2d  = CLHEP::twopi*radius;
    double arc3d  = arc2d/sz;
    double zarc   = arc2d*std::abs(ct);
    double beta   = sim.startMomentum().getV().mag()/sim.startMomentum().e();

    double t0 = sim.startGlobalTime();
    cout << "beta: " 
         << beta << " " 
         << sim.startMomentum().getV().mag() << " "
         << sim.startMomentum().e() << " "
         << radius << " "
         << zarc
         << endl;

    cout << "Constants: "
         << Constants::c <<  " " 
         << BField::cmTeslaToGeVc
         << endl;

    for ( int i=0; i<points->size(); ++i){
      
      // Aliases (references), used for readability.
      StepPointMC const& hit = (*points)[i];
      CLHEP::Hep3Vector  const& pos = hit.position();
      CLHEP::Hep3Vector  const& mom = hit.momentum();
      
      // Get the straw information, also by reference.
      Straw const&      straw = tracker.getStraw(hit.strawIndex());
      CLHEP::Hep3Vector const& mid   = straw.getMidPoint();
      CLHEP::Hep3Vector const& w     = straw.getDirection();

      // Compute straight line approximation of the drift distance.
      TwoLinePCA pca( mid, w, pos, mom);
      double dcaTrue = pca.dca();
      double dca  = dcaTrue + CLHEP::RandGaussQ::shoot(0.,sigma);

      // Arc length from poca to step point.
      double arc = (hit.time()-t0)*beta*CLHEP::c_light;
      double  turns = arc/arc3d;

      double dz     = (pos.z()-startPar(kz0));
      double turns2 = dz/radius/ct/CLHEP::twopi;

      /*
      // To test computation of poca.
      HepSymMatrix dummy(5,1);
      HelixTraj    trktraj(startPar,dummy);
      HepPoint     wirestart( mid.x(), mid.y(), mid.z());
      TrkLineTraj  wiretraj(wirestart,w);
      double       arcestimate(arc);
      double       wireestimate(straw.getHalfLength());
      TrkPoca      wpoca(trktraj,arcestimate,wiretraj,wireestimate);
      */

      // Fill diagnostic histograms and make diagnostic printout.
      if ( _diagLevel > 0 ) {

        // Check the radius of the reference point in the local
        // coordinates of the straw.  It should be 2.5 mm.
        double s = w.dot(pos-mid);
        CLHEP::Hep3Vector point = pos - (mid + s*w);
        _hCheckPointRadius->Fill(point.mag());
        _hDriftDist->Fill(pca.dca());

        if ( ncalls < _maxFullPrint ) {
          cerr << "Readback hit: "
               << event.id().event() << " "
               << i                  <<  " "
               << hit.trackId()      << "   "
               << hit.volumeId()     << " "
               << straw.Id()         << " | "
               << pca.dca()          << " "
               << pos                << " "
               << mom                << " "
               << point.mag()        << " "
               << s                  << " "
               << hit.time()-t0      << " " 
               << arc << " " 
               << turns << " " 
               << turns2 << " " 
               << endl;
        }
      }

    } // end loop over StepPointMC's

  } // end of ::analyze.

} // end of namespace mu2e

using mu2e::KalmanT01;
DEFINE_FWK_MODULE(KalmanT01);
