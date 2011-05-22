//
// Module to understand how to use the BaBar Kalman filter package.
// Not for general use.
//
// $Id: KalmanT01_module.cc,v 1.6 2011/05/22 20:28:13 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/22 20:28:13 $
//
// Original author Rob Kutschke
//
// Notes:
// 1) First version is designed to work with measurements only; ie
//    no multiple scattering.  It should work with the
//

// C++ includes.
#include <cmath>
#include <iostream>
#include <string>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Root includes.
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TNtuple.h"

// Mu2e includes.
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "KalmanTests/inc/printTrkParams.hh"
#include "LTrackerGeom/inc/LTracker.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
#include "Mu2eUtilities/inc/toHepPoint.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "ToyDP/inc/GenParticleCollection.hh"
#include "ToyDP/inc/SimParticleCollection.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"

// Babar Kalman filter includes
#include "BField/BFieldFixed.hh"
#include "TrkBase/TrkHelixUtils.hh"

// Other includes.
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Random/RandGaussQ.h"

using namespace std;

namespace mu2e {

  // Indices of track parametes within the
  // For the () operator;
  // The [] operator uses 0 based numbers.  See TrkBase/HelixTraj.hh
  // Is there an official Babar version of this for 1 based indices?
  enum TrkParIdx { kd0 =1, kphi0, kom, kz0, kct };

  class KalmanT01 : public art::EDAnalyzer {
  public:
    explicit KalmanT01(fhicl::ParameterSet const& pset) :

      // Parameters
      _generatorModuleLabel(pset.get<string>("generatorModuleLabel","generate")),
      _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel", "g4run")),
      _diagLevel(pset.get<int>("diagLevel",1)),
      _maxFullPrint(pset.get<int>("maxFullPrint",5)),
      _hDriftDist(0),
      _hCheckPointRadius(0),

      // Random number distributions
      _gaussian( createEngine( get_seed_value(pset)) ),

      _messageCategory("FitTest"){
    }
    virtual ~KalmanT01() { }

    virtual void beginJob();

    void analyze( const art::Event& e);

  private:

    // Name of the module that made these hits.
    string _generatorModuleLabel;

    // Diagnostics level.
    int _diagLevel;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    // Random number distributions.
    CLHEP::RandGaussQ _gaussian;

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

  void KalmanT01::beginJob(){

    // Create histograms if diagnostics are enabled.
    if ( _diagLevel > 0 ){

      art::ServiceHandle<art::TFileService> tfs;

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

  void KalmanT01::analyze(const art::Event& event) {

    // Counter used by debug printout.
    static int ncalls(0);
    ++ncalls;

    // Tracker origin in the Mu2e system.
    static CLHEP::Hep3Vector trackerOrigin( -3904., 0, 10200.);

    // Get a reference to one of the L or T trackers.  Throw exception if not successful.
    const Tracker& tracker = getTrackerOrThrow();

    // Ask the event to give us a handle to the requested hits.
    art::Handle<StepPointMCCollection> points;
    static const string collectionName("tracker");
    event.getByLabel(_generatorModuleLabel,collectionName,points);

    // Get handles to the generated and simulated particles.
    art::Handle<GenParticleCollection> gensHandle;
    event.getByLabel(_generatorModuleLabel, gensHandle);
    const GenParticleCollection& gens = *gensHandle;

    art::Handle<SimParticleCollection> simsHandle;
    event.getByLabel(_g4ModuleLabel, simsHandle);
    const SimParticleCollection& sims = *simsHandle;

    // Eventually get this from a calibration DB.
    static const double sigma=0.1*CLHEP::mm;

    if ( sims.empty() ) {
      cerr << "No hits. Skipping event." << endl;
      return;
    }

    const SimParticle& sim = sims.begin()->second;

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
      double dca  = dcaTrue + _gaussian.fire(0.,sigma);

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
               << straw.id()         << " | "
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
DEFINE_ART_MODULE(KalmanT01);
