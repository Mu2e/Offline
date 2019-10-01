//
// An EDProducer Module that checks radiative pi decays
//
// $Id: RPC_module.cc,v 1.15 2013/10/21 20:44:04 genser Exp $
// $Author: genser $
// $Date: 2013/10/21 20:44:04 $
//
// Original author R. Bernstein
//

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/Randomize.h"
#include "GeneralUtilities/inc/RootNameTitleHelper.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TMLPAnalyzer.h"
#include "TMath.h"
#include "TMultiLayerPerceptron.h"
#include "TNtuple.h"
#include "TSpectrum.h"
#include "TSpectrum2.h"
#include "TSpectrum3.h"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "cetlib_except/exception.h"
#include "cetlib/pow.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <cmath>
#include <ctime>
#include <iostream>
#include <set>
#include <string>
#include <utility>

using CLHEP::Hep3Vector;
using CLHEP::RandPoisson;
using cet::diff_of_squares;

using namespace mu2e;
using namespace std;

namespace mu2e {

  //--------------------------------------------------------------------
  //
  //

  class RPC : public art::EDAnalyzer{
  public:

    typedef SimParticleCollection::key_type key_type;

    RPC(fhicl::ParameterSet const& pset):
      art::EDAnalyzer(pset),

      //
      // Run time parameters
      _g4ModuleLabel(pset.get<string>("g4ModuleLabel")),
      _trackerStepPoints(pset.get<string>("trackerStepPoints","tracker")),
      _minimumEnergy(pset.get<double>("minimumEnergy")),
      _maxFullPrint(pset.get<int>("maxFullPrint",10)),
      _nAnalyzed(0),
      _messageCategory("RPC"),
      _dEdXelow(0.),
      _dEdXehi(20.),
      _dEdXnbins(2000) {}
    virtual ~RPC() {}

    virtual void beginJob();
    virtual void endJob();

    virtual void beginRun(art::Run const &r );

    virtual void beginSubRun(art::SubRun const& lblock );

    // This is called for each event.
    void analyze(const art::Event& e );


  private:

    // Module label of the g4 module that made the hits.
    std::string _g4ModuleLabel;

    // Name of the tracker StepPoint collection
    std::string _trackerStepPoints;

    // Cut on the minimum energy.
    double _minimumEnergy;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;


    // Number of events analyzed.
    int _nAnalyzed;


    // A category for the error logger.
    const std::string _messageCategory;


    TH1D* _piCaptureConvertedElectronMomentum;
    TH1D* _piCaptureConvertedElectronMomentumSignal;
    TH1D* _piCaptureConvertedElectronCosTheta;
    TH1D* _piCaptureConvertedElectronCosThetaSignal;

    TH1D* _piCaptureConvertedElectronMomentumWeighted;
    TH1D* _piCaptureConvertedElectronMomentumWeightedSignal;
    TH1D* _piCaptureConvertedElectronCosThetaWeighted;
    TH1D* _piCaptureConvertedElectronCosThetaWeightedSignal;

    TH1D* _piCaptureConvertedElectronMomentumHitTracker;
    TH1D* _piCaptureConvertedElectronMomentumSignalHitTracker;
    TH1D* _piCaptureConvertedElectronMomentumAtTracker;
    TH1D* _piCaptureConvertedElectronMomentumSignalAtTracker;
    TH1D* _piCaptureConvertedElectronCosThetaHitTracker;
    TH1D* _piCaptureConvertedElectronCosThetaSignalHitTracker;
    TH1D* _piCaptureConvertedElectronCosThetaAtTracker;
    TH1D* _piCaptureConvertedElectronCosThetaSignalAtTracker;
    TH1D* _piCaptureMomentumEnteringTracker;
    TH1D* _piCaptureMomentumLostBeforeTracker;


    TH1D* _piCaptureConvertedPositronMomentum;
    TH1D* _piCaptureConvertedPositronMomentumSignal;
    TH1D* _piCaptureConvertedPositronCosTheta;
    TH1D* _piCaptureConvertedPositronCosThetaSignal;

    TH1D* _piCaptureConvertedFinalSpectrum;
    TH1D* _piCaptureConvertedFinalSpectrumWeighted;

    TH1D* _conversionAsymmetry;
    TH1D* _numberOfHitStraws;
    TH1D* _piCaptureZofHit;
    TH1D* _piCaptureZDiff;

    //
    // next block handles energy loss in stopping foils and proton absorber
    //
    // dE/dx spectrum as a continuous function.
    const double energyLossSpectrum(const double e);
    const double isotropicEnergyLossSpectrum(const double e);
    // Compute a binned representation of the dE/dx spectrm.
    std::vector<double> binnedEnergyLossSpectrum();
    std::vector<double> isotropicBinnedEnergyLossSpectrum();
    double _dEdXelow;         //< lower dE/dx for binned plot
    double _dEdXehi;          //< upper
    int    _dEdXnbins;        //< number of bins
    TH1D* _piCaptureEnergyLossSpectrum;
    //CLHEP::RandGeneral _dEdXspectrum;


    void bookEventHistos(double const elow, double const ehigh);
    void fillEventHistos();

  };


  void RPC::beginJob( ){

    // Get access to the TFile service.
    //    art::ServiceHandle<art::TFileService> tfs;

  }

  void RPC::endJob(){


    cout << " time for this job was:  " << clock()/CLOCKS_PER_SEC << endl;

  }



  void RPC::beginRun(art::Run const& run ){
  }

  void RPC::beginSubRun(art::SubRun const& lblock ){
  }


  void RPC::analyze(const art::Event& event ) {

    static int ncalls(0);
    ++ncalls;

    // Maintain a counter for number of events seen.
    ++_nAnalyzed;

    //    if (ncalls >= 81000){cout << "ncalls = " << ncalls << endl;} //assert(2==1);


    // Book histogram on the first call regardless
    double elow = 103.5;
    double ehigh = 105.;

    if ( ncalls == 1){
      bookEventHistos(elow,ehigh);
      //
      // fill energy loss histo for future use
      vector<double> energyLossVector = isotropicBinnedEnergyLossSpectrum();
      double binsize = (_dEdXehi - _dEdXelow)/_dEdXnbins;
      for (int i = 0; i<=_dEdXnbins; ++i){
        _piCaptureEnergyLossSpectrum->Fill(_dEdXelow + i*binsize + binsize/2.,energyLossVector[i]);
      }
      //
      // print integral as check
      double energyLossIntegral = _piCaptureEnergyLossSpectrum->Integral();
      cout<< "integral of energy loss spectrum = " << energyLossIntegral << endl;
      _piCaptureEnergyLossSpectrum->Scale(1./energyLossIntegral);

    }
    mf::LogVerbatim log(_messageCategory);
    log << "RPC event #: "
        << event.id();

    //
    // start looking through SimParticles
    art::Handle<SimParticleCollection> simParticles;
    event.getByLabel(_g4ModuleLabel, simParticles);

    // Handle to information about G4 physical volumes.
    art::Handle<PhysicalVolumeInfoCollection> volumes;
    event.getRun().getByLabel(_g4ModuleLabel, volumes);

    // Some files might not have the SimParticle and volume information.
    bool haveSimPart = ( simParticles.isValid() && volumes.isValid() );

    // Other files might have empty collections.
    if ( haveSimPart ){
      haveSimPart = !(simParticles->empty() || volumes->empty());
    }

    //
    // get handle to hits:  with _trackerStepPoints, if there are hits they belong to the tracker
    art::Handle<StepPointMCCollection> hits;
    event.getByLabel(_g4ModuleLabel,_trackerStepPoints,hits);


    double photonEnergy = 0.;
    double electronEnergy = 0.;
    double positronEnergy = 0.;

    CLHEP::Hep3Vector  momentumAtEntranceToTracker = CLHEP::Hep3Vector();

    if (haveSimPart){
      for ( SimParticleCollection::const_iterator ithPart=simParticles->begin();
            ithPart!=simParticles->end(); ++ithPart ){

        SimParticle const& sim = ithPart->second;
        PhysicalVolumeInfo const& startVol = volumes->at(sim.startVolumeIndex());
        //
        // is this the initial photon?
        if ( !sim.hasParent() ){
          if (sim.pdgId() != PDGCode::gamma) {
            //
            // this can't happen if we're studying RPCs so throw and die
            throw cet::exception("GEOM")
              << "RPC with a parent not a photon, but parent PDG code is "
              << sim.pdgId();
          } else
            {   CLHEP::HepLorentzVector photonMomentum = sim.startMomentum();
              photonEnergy = photonMomentum.e();
            }

        }
        //      cout << " volumename = " << startVol.name() << endl;
        //
        // check three things:  (1) the mother is the original photon, (2) you're an e+ or e-, and (3) the photon converts in the foil
        if ( sim.parentId().asInt() == 1 && startVol.name().compare(0,11,"TargetFoil_") == 0 ){

          if (sim.pdgId() == PDGCode::e_minus) {
            bool electronHitTracker = false;
            bool electronAccepted = false;
            //bool hitEnoughStraws = false;
            double zmin = +99999.;
            double zmax = -99999.;

            const CLHEP::HepLorentzVector& electronMomentum = sim.startMomentum();
            electronEnergy = electronMomentum.e();

            //
            // i want the momentum at the tracker, not the momentum at birth.  RPCs can smear
            // in and out of the signal window.  loop over the hits and find the right one
            bool firstHitOnElectronTrack = true;
            for( size_t i=0; i<hits->size(); ++i ){

              const StepPointMC& hit = (*hits)[i];

              //step point mc associated with generated electrons
              key_type trackId = hit.trackId();
              //
              // now I have the track Id of this hit.
              SimParticle const& simParticleHit = simParticles->at(trackId);
              //
              // is this simParticle associated with the hit the same as the one I started with?
              if (sim.id() == simParticleHit.id()){
                electronHitTracker = true;

                //
                //let's make a distribution of the z position of all hits associated with this track
                _piCaptureZofHit->Fill(hit.position().z());

                //
                // and my trigger will be that there's 1.5 meters between first and last z-hit on this track.
                if (zmin > hit.position().z()) {zmin = hit.position().z();}
                if (zmax < hit.position().z()) {zmax = hit.position().z();}
                //          cout << "z of hit " << hit.position().z() << endl;

                //
                // the momentum at the first hit of the tracker is the momentum entering the tracker
                if (firstHitOnElectronTrack){
                  firstHitOnElectronTrack = false;
                  momentumAtEntranceToTracker = hit.momentum();
                  //double momentum = sqrt(diff_of_squares(electronMomentum.e(), electronMomentum.invariantMass()));
                  //                  cout << "momentum at entrance to tracker = " << momentumAtEntranceToTracker.mag() << endl;
                  //cout << "original momentum was           = " << momentum << endl;
                  //cout << "track Id                        = " << trackId  << endl;
                  if (momentumAtEntranceToTracker.mag() < 1.0) {cout << " bad one!!" << endl;}
                }

              }

            } // end loop over hits.
            double zDiffCut = 1750.; // very arbitrary choice and bad hardwired numbers.  see plot of zdiff for rationale.  cutting short tracks and ones that enter toward back.  used conversion electrons to make this choice
            if (electronHitTracker){
              _piCaptureZDiff->Fill(zmax - zmin);
              //              cout << "zmin,zmax " << zmin << " " << zmax << endl;
              //
              // key line that defines electron acceptance
              if ( (zmax - zmin) > zDiffCut && abs(zmax - zmin) < 4000.) {
                electronAccepted = true;
              }
            }
            //
            // the above doesn't give a great angular distribution;
            //but just hitting the tracker looks good up to an overall scale
            //
            // uncomment iff i'm forcing the normalization to 19% from TTracker memo anyway, and using MECO's weight function
            //electronAccepted = false;
            //if (electronHitTracker) {electronAccepted = true;}
            //
            // can't demand particle enters tracker if I apply the weight function below.  But I do want to know what
            // the energy loss distribution is for electrons that do hit.  So:
            double momentum = sqrt(diff_of_squares(electronMomentum.e(), electronMomentum.invariantMass()));
            if (!electronAccepted){
              momentumAtEntranceToTracker = electronMomentum.getV();
              //              double zcheck = electronMomentum.getZ();
              //cout << "checking .vect " << momentumAtEntranceToTracker << "\n" << electronMomentum << endl;
            }

            double cost = electronMomentum.cosTheta();

            //
            // all electrons, broad spectrum, no weights or cuts
            _piCaptureConvertedElectronMomentum->Fill(momentum);
            _piCaptureConvertedElectronCosTheta->Fill(cost);
            if (momentum > elow && momentum < ehigh){
              _piCaptureConvertedElectronMomentumSignal->Fill(momentum);
              _piCaptureConvertedElectronCosThetaSignal->Fill(cost);
            }

            //
            // here is a weight function from a fit, documented separately.
            // the idea is I have an acceptance function from GMC, and then weight the electrons by their acceptance
            // as a function of cos theta
            double weight = 0.;

            //
            // the integrals of these weight functions over d(cos theta) from 0 to 1 and -1 to 0 are .150 and .133 respectively,
            // see doc-db 1087. weight integral over -1 to 1 to 19% to match Rashid's memo on T-tracker.  Used to use 29% but that
            // was before calorimeter acceptance and track quality cuts.
            //
            // last fudge factor is because Rashid's 19% is for events in signal box.  But that includes energy loss, and I'm too early in the code
            // for that, so scale up using energy loss

            //one has to be very careful with these weights.  The weight is for 105 MeV conversion electrons, so if I look in a restricted region
            // where the acceptance is about the same as conversion electrons, it's fine.

            //
            // last factor of two is because \int_1^1 d(cos theta) = 2 so I need this to give me 19% acceptance.  I'm essentially multiplying
            // the PDF by the bin width.  This ends up being the number accepted/number generated as a function of cos theta, and that's what we want
            // for a weight.
            if (cost >= 0){
              weight = 2.*(19./29.)*(1/.04)*(1.48E-02)*TMath::Exp( - TMath::Power((cost - 0.200),2)/(2.*.19*.19));
            }
            if (cost < 0 && cost >= -0.54) {
              weight = 2.*(19./29.)*(1/.04)*( (1.04e-02) - 3.19e-02*cost - 9.3e-02*cost*cost);
           }
            if (cost < -0.54 || weight < 0.) {weight = 0.;}


            //
            // and final weight fudge, based on energy loss for
            // isotropic conversion electrons entering tracker. That will make acceptance match Rashid's memo.
            // RPCs lose more energy but too bad for them, that's not the factor I want.   Not perfect since entering tracker isn't necessarily
            // a good representation of reconstructed events, but not wrong either.  Really 1.50 within errors
            weight *= 1.5;

            //
            // all electrons but weighted
            _piCaptureConvertedElectronMomentumWeighted->Fill(momentum,weight);
            _piCaptureConvertedElectronCosThetaWeighted->Fill(cost,weight );
            if (momentum > elow && momentum < ehigh){
              _piCaptureConvertedElectronMomentumWeightedSignal->Fill(momentum,weight);
              _piCaptureConvertedElectronCosThetaWeightedSignal->Fill(cost,weight );
            }
            //
            // weight events according to energy loss spectrum (big clause below "demand they entered tracker" computes this,
            //     we re-use here
            //   there's some correlation between the weight and the angle so this isn't
            //  quite right

            //
            // smear energy loss downward according to the function simulated above in an earlier job

            double momentumForPlot = momentum-_piCaptureEnergyLossSpectrum->GetRandom();
            _piCaptureConvertedFinalSpectrumWeighted->Fill(momentumForPlot,weight);
            //
            // and now demand they entered tracker -- these should not be weighted since the weight mocks up the acceptance
            // however, this allows us to plot what we measure -- momentum at the tracker after energy loss in
            // stopping target and foils -- hence we get RPCs that are smeared down and lose RPCs that are smeared out.
            //  If the # of events and distributions in cos theta and momentum (at birth)
            // are the same as when weighted, then this hitting the tracker is a good mock-up of accepted events.
            if(electronAccepted){
              double momentumAtTracker = momentumAtEntranceToTracker.mag();
              _piCaptureConvertedFinalSpectrum->Fill(momentumAtTracker);
              //
              // this plot only has meaning if you entered the tracker
              _piCaptureMomentumLostBeforeTracker->Fill(momentum - momentumAtEntranceToTracker.mag());

              //cout << " inside plot, momentumAtTracker = " << momentumAtTracker << endl;
              _piCaptureConvertedElectronMomentumHitTracker->Fill(momentum);
              _piCaptureConvertedElectronCosThetaHitTracker->Fill(cost);
              _piCaptureConvertedElectronCosThetaAtTracker->Fill( momentumAtEntranceToTracker.cosTheta());
              _piCaptureConvertedElectronMomentumAtTracker->Fill(momentumAtTracker);

              if (momentumAtTracker > elow && momentumAtTracker < ehigh){
                _piCaptureConvertedElectronMomentumSignalHitTracker->Fill(momentum);
                _piCaptureConvertedElectronMomentumSignalAtTracker->Fill(momentumAtTracker);
                _piCaptureConvertedElectronCosThetaSignalHitTracker->Fill(cost);
                _piCaptureConvertedElectronCosThetaSignalAtTracker->Fill(momentumAtEntranceToTracker.cosTheta());
              }
            }
          }

          if (sim.pdgId() == PDGCode::e_plus) {
            CLHEP::HepLorentzVector electronMomentum = sim.startMomentum();
            positronEnergy = electronMomentum.e();
            double momentum = sqrt(diff_of_squares(electronMomentum.e(), electronMomentum.invariantMass()));
            _piCaptureConvertedPositronMomentum->Fill(momentum);
            if(momentum > elow && momentum < ehigh) _piCaptureConvertedPositronMomentumSignal->Fill(momentum);
            _piCaptureConvertedPositronCosTheta->Fill(electronMomentum.cosTheta() );
            if (momentum > elow && momentum < ehigh) _piCaptureConvertedPositronCosThetaSignal->Fill(electronMomentum.cosTheta() );
          }
        }
      }
      //
      // check geant4's photon conversion.  if there was no conversion you won't see anything so check
      if (simParticles->size() >= 3 && electronEnergy > 0. && positronEnergy > 0. && photonEnergy > 0.){
        _conversionAsymmetry->Fill( abs(electronEnergy - positronEnergy)/ (electronEnergy + positronEnergy) );
      }
    }
  }
  void RPC::bookEventHistos(double const elow, double const ehigh)
  {
    //    cout << "booking histos" << endl; assert(2==1);
    art::ServiceHandle<art::TFileService> tfs;



    _piCaptureEnergyLossSpectrum =
      tfs->make<TH1D>( "piCaptureEnergyLossSpectrum",
                       "Pi Capture Energy Loss Spectrum",_dEdXnbins,_dEdXelow,_dEdXehi);
    _piCaptureZofHit =
      tfs->make<TH1D>( "piCaptureZofHit",
                       "Pi Capture Z of Hit", 200, -2000., 2000.);


    //
    // all events, unweighted
    _piCaptureConvertedElectronMomentum =
      tfs->make<TH1D>( "piCaptureConvertedElectronMomentum",
                       "Pi Capture Converted Electron Momentum", 200, 0., 200.);
    _piCaptureConvertedFinalSpectrum =
      tfs->make<TH1D>( "piCaptureConvertedFinalSpectrum",
                       "Pi Capture Converted Final Spectrum", 15, elow, ehigh);
    _piCaptureConvertedFinalSpectrumWeighted =
      tfs->make<TH1D>( "piCaptureConvertedFinalSpectrumWeighted",
                       "Pi Capture Converted Final Spectrum, Weighted for Acceptance", 15, elow, ehigh);
    _piCaptureConvertedElectronMomentumSignal =
      tfs->make<TH1D>( "piCaptureConvertedElectronMomentumSignal",
                       "Pi Capture Converted Electron MomentumSignal", 15, elow, ehigh);
    _piCaptureConvertedElectronCosTheta =
      tfs->make<TH1D>( "piCaptureConvertedElectronCosTheta",
                       "Pi Capture Converted Electron CosTheta", 200, -1., 1.);
    _piCaptureConvertedElectronCosThetaSignal =
      tfs->make<TH1D>( "piCaptureConvertedElectronCosThetaSignal",
                       "Pi Capture Converted Electron CosThetaSignal", 200, -1., 1.);

    //
    // same, but now weighted by acceptance
    _piCaptureConvertedElectronMomentumWeighted =
      tfs->make<TH1D>( "piCaptureConvertedElectronMomentumWeighted",
                       "Pi Capture Converted Electron MomentumWeighted", 200, 0., 200.);
    _piCaptureConvertedElectronMomentumWeightedSignal =
      tfs->make<TH1D>( "piCaptureConvertedElectronMomentumWeightedSignal",
                       "Pi Capture Converted Electron MomentumWeightedSignal", 15, elow, ehigh);
    _piCaptureConvertedElectronCosThetaWeighted =
      tfs->make<TH1D>( "piCaptureConvertedElectronCosThetaWeighted",
                       "Pi Capture Converted Electron CosThetaWeighted", 200, -1., 1.);
    _piCaptureConvertedElectronCosThetaWeightedSignal =
      tfs->make<TH1D>( "piCaptureConvertedElectronCosThetaWeightedSignal",
                       "Pi Capture Converted Electron CosThetaWeightedSignal", 200, -1., 1.);


    //
    // and now demand a tracker hit -- but not weighted, since tracker hit is like an acceptance is like a weight
    _piCaptureConvertedElectronMomentumHitTracker =
      tfs->make<TH1D>( "piCaptureConvertedElectronMomentumHitTracker",
                       "Pi Capture Converted Electron Momentum and Hit Tracker", 200, 0., 200.);
    //
    // this plot is the original momentum of electrons that end up in signal region after energy loss so ehigh
    // needs to be bigger
    _piCaptureConvertedElectronMomentumSignalHitTracker =
      tfs->make<TH1D>( "piCaptureConvertedElectronMomentumSignalHitTracker",
                       "Pi Capture Converted Electron Momentum Signal and Hit Tracker", 15, elow, 110.);
    _piCaptureConvertedElectronCosThetaHitTracker =
      tfs->make<TH1D>( "piCaptureConvertedElectronCosThetaHitTracker",
                       "Pi Capture Converted Electron CosTheta Hit Tracker", 200, -1., 1.);
    _piCaptureConvertedElectronCosThetaSignalHitTracker =
      tfs->make<TH1D>( "piCaptureConvertedElectronCosThetaSignalHitTracker",
                       "Pi Capture Converted Electron CosTheta Signal Hit Tracker", 200, -1., 1.);
    _piCaptureConvertedElectronCosThetaAtTracker=
      tfs->make<TH1D>( "piCaptureConvertedElectronCosThetaAtTracker",
                       "Pi Capture Converted Electron CosTheta At Tracker", 200, -1., 1.);
    _piCaptureConvertedElectronCosThetaSignalAtTracker =
      tfs->make<TH1D>( "_piCaptureConvertedElectronCosThetaSignalAtTracker",
                       "Pi Capture Converted Electron CosTheta Signal At Tracker", 200, -1., 1.);
    _piCaptureConvertedElectronMomentumAtTracker =
      tfs->make<TH1D>("piCaptureConvertedElectronMomentumAtTracker","Electron Momentum At Tracker", 200,0.,200.);
    _piCaptureConvertedElectronMomentumSignalAtTracker =
      tfs->make<TH1D>("piCaptureConvertedElectronMomentumSignalAtTracker","Electron Momentum Signal At Tracker", 15,elow,ehigh);
    _piCaptureMomentumLostBeforeTracker =
      tfs->make<TH1D>("_piCaptureMomentumLostBeforeTracker","Electron Momentum Lost Before Tracker", 100,0.,10.);



    //
    //positron plots
    _piCaptureConvertedPositronCosTheta =
      tfs->make<TH1D>( "piCaptureConvertedPositronCosTheta",
                       "Pi Capture Converted Positron CosTheta", 200, -1., 1.);
    _piCaptureConvertedPositronCosThetaSignal =
      tfs->make<TH1D>( "piCaptureConvertedPositronCosThetaSignal",
                       "Pi Capture Converted Positron CosThetaSignal", 200, -1., 1.);
    _piCaptureConvertedPositronMomentum =
      tfs->make<TH1D>( "piCaptureConvertedPositronMomentum",
                       "Pi Capture Converted Positron Momentum", 200, 0., 200.);
    _piCaptureConvertedPositronMomentumSignal =
      tfs->make<TH1D>( "piCaptureConvertedPositronMomentumSignal",
                       "Pi Capture Converted Positron MomentumSignal", 20, elow, ehigh);



    _conversionAsymmetry =
      tfs->make<TH1D>("conversionAsymmetry","electron-positron/ephoton", 100,0.,1.);

    _numberOfHitStraws =
      tfs->make<TH1D>("numberOfHitStraws"," Number Of Hit Straws", 100,0.,100.);

    _piCaptureZDiff =
      tfs->make<TH1D>("piCaptureZDiff","Length of Track", 100,0.,5000.);


    //
    // set up dE/dx spectrum

  } //bookEventHistos


  const double RPC::energyLossSpectrum(const double e){
    //
    // all this from fitting energy loss plotted elsewhere in this job and fitting
    double loss = 0.;
    if (e <= 2.5)
      { loss = 2177.6*TMath::Landau(e,1.25,.37463);}// fudge numbers from histo in denom
    if (e >= 2.5 && e <= 10) {loss = TMath::Exp(6.40 + -0.5899*(e) + .02638*e*e); }
    if (e > 10) {loss = 0.;}
    if (e<0){
      throw cet::exception("RANGE")
        << "Nonsense energy in RPC_plugin.cc="
        << e
        << "\n";
    }


    // normalize loss histo to unity (above depended on statistics of particular job)
    // this is done anyway in GetRandom method, just useful to do it here
    return loss/10032.;
  }

  // Compute a binned representation of the energy loss spectrum.
  std::vector<double> RPC::binnedEnergyLossSpectrum(){

    // Sanity check.
    if (_dEdXnbins <= 0) {
      throw cet::exception("RANGE")
        << "Nonsense RPC_plugin.nbins requested="
        << _dEdXnbins
        << "\n";
    }

    // Bin width.
    double dE = (_dEdXehi - _dEdXelow) / _dEdXnbins;

    // Vector to hold the binned representation of the energy spectrum.
    std::vector<double> spectrum;
    spectrum.reserve(_dEdXnbins);
    double sum = 0.;
    for (int ib=0; ib<_dEdXnbins; ib++) {
      double x = _dEdXelow+(ib+0.5) * dE;
      sum += x;
      spectrum.push_back(energyLossSpectrum(x));
    }
    cout << "summed energy loss " << sum << endl;
    return spectrum;
  } // RPC::binnedEnergyLossSpectrum



  const double RPC::isotropicEnergyLossSpectrum(const double e){
    //
    // all this from fitting energy loss plotted elsewhere in this job and fitting
    double loss = 0.;
    if (e <= 4.)
      { loss = 2567.38*TMath::Landau(e,0.89539,.202878);}// fudge numbers from histo in denom
    if (e > 4.0 && e <= 10) {loss = TMath::Exp(4.5029+ -0.6358*(e) + 2.619E-02*(e)*(e)); }
    if (e>10){loss = 0.;}
    if (e<0){
      throw cet::exception("RANGE")
        << "Nonsense energy in RPC_plugin.cc="
        << e
        << "\n";
    }


    // normalize loss histo to unity (above depended on statistics of particular job)
    // done later anyway in GetRandom() but nice to have here
    return loss/5170.;
  }

  // Compute a binned representation of the energy loss spectrum.
  std::vector<double> RPC::isotropicBinnedEnergyLossSpectrum(){

    // Sanity check.
    if (_dEdXnbins <= 0) {
      throw cet::exception("RANGE")
        << "Nonsense RPC_plugin.nbins requested="
        << _dEdXnbins
        << "\n";
    }

    // Bin width.
    double dE = (_dEdXehi - _dEdXelow) / _dEdXnbins;

    // Vector to hold the binned representation of the energy spectrum.
    std::vector<double> spectrum;
    spectrum.reserve(_dEdXnbins);
    double sum = 0.;
    for (int ib=0; ib<_dEdXnbins; ib++) {
      double x = _dEdXelow+(ib+0.5) * dE;
      sum += x;
      spectrum.push_back(isotropicEnergyLossSpectrum(x));
    }
    cout << "summed energy loss " << sum << endl;
    return spectrum;
  } // RPC::binnedEnergyLossSpectrum


} // namespace mu2e


DEFINE_ART_MODULE(RPC);
