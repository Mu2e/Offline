//
// An EDProducer Module that checks conversion electrons
//
// $Id: CEL_module.cc,v 1.15 2013/10/21 20:44:04 genser Exp $
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
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "cetlib_except/exception.h"
#include "cetlib/pow.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <cmath>
#include <cstdlib>
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
  class CEL;
}

//--------------------------------------------------------------------
//
//

class mu2e::CEL : public art::EDAnalyzer{
public:

  typedef SimParticleCollection::key_type key_type;

  CEL(fhicl::ParameterSet const& pset):
    art::EDAnalyzer(pset),
    //
    // Run time parameters
    _g4ModuleLabel(pset.get<string>("g4ModuleLabel")),
    _trackerStepPoints(pset.get<string>("trackerStepPoints","tracker")),
    _minimumEnergy(pset.get<double>("minimumEnergy")),
    _maxFullPrint(pset.get<int>("maxFullPrint",10)),
    _nAnalyzed(0),
    _messageCategory("CEL"),
    _dEdXelow(0.),
    _dEdXehi(20.),
    _dEdXnbins(2000) {}
  virtual ~CEL() {}

  virtual void beginJob();
  virtual void endJob();

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


  TH1D* _cELConvertedElectronMomentum;
  TH1D* _cELConvertedElectronMomentumSignal;
  TH1D* _cELConvertedElectronCosTheta;
  TH1D* _cELConvertedElectronCosThetaSignal;

  TH1D* _cELConvertedElectronMomentumWeighted;
  TH1D* _cELConvertedElectronMomentumWeightedSignal;
  TH1D* _cELConvertedElectronCosThetaWeighted;
  TH1D* _cELConvertedElectronCosThetaWeightedSignal;

  TH1D* _cELConvertedElectronMomentumHitTracker;
  TH1D* _cELConvertedElectronMomentumSignalHitTracker;
  TH1D* _cELConvertedElectronMomentumAtTracker;
  TH1D* _cELConvertedElectronMomentumSignalAtTracker;
  TH1D* _cELConvertedElectronCosThetaHitTracker;
  TH1D* _cELConvertedElectronCosThetaSignalHitTracker;
  TH1D* _cELConvertedElectronCosThetaAtTracker;
  TH1D* _cELConvertedElectronCosThetaSignalAtTracker;
  TH1D* _cELMomentumEnteringTracker;
  TH1D* _cELMomentumLostBeforeTracker;


  TH1D* _cELConvertedPositronMomentum;
  TH1D* _cELConvertedPositronMomentumSignal;
  TH1D* _cELConvertedPositronCosTheta;
  TH1D* _cELConvertedPositronCosThetaSignal;

  TH1D* _cELConvertedFinalSpectrum;
  TH1D* _cELConvertedFinalSpectrumWeighted;

  TH1D* _conversionAsymmetry;
  TH1D* _numberOfHitStraws;
  TH1D* _cELZofHit;
  TH1D* _cELZDiff;

  //
  // next block handles energy loss in stopping foils and proton absorber
  //
  // dE/dx spectrum as a continuous function.
  const double energyLossSpectrum(const double e);
  // Compute a binned representation of the dE/dx spectrm.
  std::vector<double> binnedEnergyLossSpectrum();
  double _dEdXelow;         //< lower dE/dx for binned plot
  double _dEdXehi;          //< upper
  int    _dEdXnbins;        //< number of bins
  TH1D* _cELEnergyLossSpectrum;
  //CLHEP::RandGeneral _dEdXspectrum;


  void bookEventHistos(double const elow, double const ehigh);
  void fillEventHistos();

};


void CEL::beginJob( ){

  // Get access to the TFile service.
  //    art::ServiceHandle<art::TFileService> tfs;

}

void CEL::endJob(){

  cout << " time for this job was:  " << clock()/CLOCKS_PER_SEC << endl;
}

void CEL::analyze(const art::Event& event ) {

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
    vector<double> energyLossVector = binnedEnergyLossSpectrum();
    double binsize = (_dEdXehi - _dEdXelow)/_dEdXnbins;
    for (int i = 0; i<=_dEdXnbins; ++i){
      _cELEnergyLossSpectrum->Fill(_dEdXelow + i*binsize + binsize/2.,energyLossVector[i]);
    }
    //
    // print integral as check
    double energyLossIntegral = _cELEnergyLossSpectrum->Integral();
    cout<< "integral of energy loss spectrum = " << energyLossIntegral << endl;
    _cELEnergyLossSpectrum->Scale(1./energyLossIntegral);

  }
  mf::LogVerbatim log(_messageCategory);
  log << "CEL event #: "
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



  CLHEP::Hep3Vector  momentumAtEntranceToTracker = CLHEP::Hep3Vector();
  cout << "have sim particle" << endl;
  if (haveSimPart){
    for ( SimParticleCollection::const_iterator i=simParticles->begin();
          i!=simParticles->end(); ++i ){

      SimParticle const& sim = i->second;

      PhysicalVolumeInfo const& startVol = volumes->at(sim.startVolumeIndex());
      //
      // is this the initial electron?
      if ( !sim.hasParent() ){
        if (sim.pdgId() != PDGCode::e_minus) {
          //
          // this can't happen if we're studying CELs so throw and die
          throw cet::exception("GEOM")
            << "CEL with a parent not an electron, but parent PDG code is "
            << sim.pdgId();
        } else
          {   CLHEP::HepLorentzVector electronMomentum = sim.startMomentum();
            double electronEnergy = electronMomentum.e();
            cout << " this had better be 105: " << electronEnergy << endl;
          }

      }
      //
      // check you're the original electron and you were born in the foil
      if ( !sim.hasParent() 
           && startVol.name().compare(0,11,"TargetFoil_") == 0 
           && sim.pdgId() == PDGCode::e_minus){
        bool electronHitTracker = false;
        bool electronAccepted = false;
        //bool hitEnoughStraws = false;
        double zmin = +99999.;
        double zmax = -99999.;

        const CLHEP::HepLorentzVector& electronMomentum = sim.startMomentum();
        //double electronEnergy = electronMomentum.e();
        //cout << "again this had better be 105. " << electronMomentum.e() << endl;

        //
        // i want the momentum at the tracker, not the momentum at birth.  CELs can smear
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
            _cELZofHit->Fill(hit.position().z());

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
        double zDiffCut = 1750.; // very arbitrary choice and bad hardwired numbers.  see plot of zdiff for rationale; cutting very short tracks
        if (electronHitTracker){
          _cELZDiff->Fill(zmax - zmin);
          //              cout << "zmin,zmax " << zmin << " " << zmax << endl;
          //
          // key line that defines electron acceptance
          if ( (zmax - zmin) > zDiffCut && abs(zmax - zmin) < 4000.) {
            electronAccepted = true;
          }
        }

        //
        // this definition of acceptance doesn't work well, so use weighting functions
        //          electronAccepted = false;
        //if (electronHitTracker){electronAccepted = true;}
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
        _cELConvertedElectronMomentum->Fill(momentum);
        _cELConvertedElectronCosTheta->Fill(cost);
        if (momentum > elow && momentum < ehigh){
          _cELConvertedElectronMomentumSignal->Fill(momentum);
          _cELConvertedElectronCosThetaSignal->Fill(cost);
        }

        //
        // here is a weight function from a fit, documented separately.
        // the idea is I have an acceptance function from GMC, and then weight the electrons by their acceptance
        // as a function of cos theta
        double weight = 0.;

        //
        // the integrals of these weight functions over d(cos theta)from 0 to 1 and from 0 to -1
        // are .150 and .133 respectively. Need to double to make integral from -1 to 1 match.
        //
        // See doc-db 1087: needs to end up with 19% in signal box so energy loss is the
        // source of the final fudge factor -- we don't know what it looks like at this point in the code.


        //
        // last factor of two is because \int_1^1 d(cos theta) = 2 so I need this to give me 19% acceptance.  I'm essentially multiplying
        // the PDF by the bin width.  This ends up being the number accepted/number generated as a function of cos theta, and that's what we want
        // for a weight

        if (cost >= 0){
          weight = 2.*(19./29.)*(1/.04)*(1.48E-02)*TMath::Exp( - TMath::Power((cost - 0.200),2)/(2.*.19*.19));
          //             cout << "cos, weight for + " << cost << " " << weight << endl;
          //if (weight < 0){cout <<"negative weight cost > 0: " << cost << " " << weight << endl;assert(2==1);}
          //weight  = 0.0;
        }
        if (cost < 0 && cost >= -0.5) {
          weight = 2.*(19./29.)*(1/.04)*( (1.04e-02) - 3.19e-02*cost - 9.3e-02*cost*cost);

          //if (weight < 0){cout <<"negative weight cost < 0: " << cost << " " << weight << endl;assert(2==1);}

          //              cout << "cos, weight for - " << cost << " " << weight << endl;
        }
        if (cost < -0.5 || weight < 0.) {weight = 0.;}

        //
        // and final weight fudge, based on energy loss for
        //isotropic conversion electrons entering tracker.  Not perfect since entering tracker isn't necessarily
        // a good representation of reconstructed events, but not wrong either.  Really 1.50 within errors
        weight *= 1.5;


        //
        // all electrons but weighted
        _cELConvertedElectronMomentumWeighted->Fill(momentum,weight);
        _cELConvertedElectronCosThetaWeighted->Fill(cost,weight );
        cout << "costheta, weight = " << cost << " " << weight << endl;
        if (momentum > elow && momentum < ehigh){
          _cELConvertedElectronMomentumWeightedSignal->Fill(momentum,weight);
          _cELConvertedElectronCosThetaWeightedSignal->Fill(cost,weight );
        }
        //
        // weight events according to energy loss spectrum (big clause below "demand they entered tracker" computes this, we re-use here
        //   there's some correlation between the weight and the angle so this isn't
        //  quite right

        double momentumForPlot = momentum-_cELEnergyLossSpectrum->GetRandom();
        _cELConvertedFinalSpectrum->Fill(momentumForPlot);
        _cELConvertedFinalSpectrumWeighted->Fill(momentumForPlot,weight);
        //
        // and now demand they entered tracker -- these should not be weighted since the weight mocks up the acceptance
        // however, this allows us to plot what we measure -- momentum at the tracker after energy loss in
        // stopping target and foils -- hence we get CELs that are smeared down and lose CELs that are smeared out.
        //  If the # of events and distributions in cos theta and momentum (at birth)
        // are the same as when weighted, then this hitting the tracker is a good mock-up of accepted events.
        if(electronAccepted){
          double momentumAtTracker = momentumAtEntranceToTracker.mag();
          //
          // this plot only has meaning if you entered the tracker
          _cELMomentumLostBeforeTracker->Fill(momentum - momentumAtEntranceToTracker.mag());
          //
          // smear energy loss downward according to the function simulated above in an earlier job


          //cout << " inside plot, momentumAtTracker = " << momentumAtTracker << endl;
          _cELConvertedElectronMomentumHitTracker->Fill(momentum);
          _cELConvertedElectronCosThetaHitTracker->Fill(cost);
          _cELConvertedElectronCosThetaAtTracker->Fill( momentumAtEntranceToTracker.cosTheta());
          _cELConvertedElectronMomentumAtTracker->Fill(momentumAtTracker);

          if (momentumAtTracker > elow && momentumAtTracker < ehigh){
            _cELConvertedElectronMomentumSignalHitTracker->Fill(momentum);
            _cELConvertedElectronMomentumSignalAtTracker->Fill(momentumAtTracker);
            _cELConvertedElectronCosThetaSignalHitTracker->Fill(cost);
            _cELConvertedElectronCosThetaSignalAtTracker->Fill(momentumAtEntranceToTracker.cosTheta());
          }
        }
      }
    }
  }
}
void CEL::bookEventHistos(double const elow, double const ehigh)
{
  //    cout << "booking histos" << endl; assert(2==1);
  art::ServiceHandle<art::TFileService> tfs;



  _cELEnergyLossSpectrum =
    tfs->make<TH1D>( "cELEnergyLossSpectrum",
                     "Conversion Electron Energy Loss Spectrum",_dEdXnbins,_dEdXelow,_dEdXehi);
  _cELZofHit =
    tfs->make<TH1D>( "cELZofHit",
                     "Conversion Electron Z of Hit", 200, -2000., 2000.);


  //
  // all events, unweighted
  _cELConvertedElectronMomentum =
    tfs->make<TH1D>( "cELConvertedElectronMomentum",
                     "Conversion Electron Converted Electron Momentum", 200, 0., 200.);
  _cELConvertedFinalSpectrum =
    tfs->make<TH1D>( "cELConvertedFinalSpectrum",
                     "Conversion Electron Converted Final Spectrum", 15, elow, ehigh);
  _cELConvertedFinalSpectrumWeighted =
    tfs->make<TH1D>( "cELConvertedFinalSpectrumWeighted",
                     "Conversion Electron Converted Final Spectrum, Weighted for Acceptance", 15, elow, ehigh);
  _cELConvertedElectronMomentumSignal =
    tfs->make<TH1D>( "cELConvertedElectronMomentumSignal",
                     "Conversion Electron Converted Electron MomentumSignal", 15, elow, ehigh);
  _cELConvertedElectronCosTheta =
    tfs->make<TH1D>( "cELConvertedElectronCosTheta",
                     "Conversion Electron Converted Electron CosTheta", 200, -1., 1.);
  _cELConvertedElectronCosThetaSignal =
    tfs->make<TH1D>( "cELConvertedElectronCosThetaSignal",
                     "Conversion Electron Converted Electron CosThetaSignal", 200, -1., 1.);

  //
  // same, but now weighted by acceptance
  _cELConvertedElectronMomentumWeighted =
    tfs->make<TH1D>( "cELConvertedElectronMomentumWeighted",
                     "Conversion Electron Converted Electron MomentumWeighted", 200, 0., 200.);
  _cELConvertedElectronMomentumWeightedSignal =
    tfs->make<TH1D>( "cELConvertedElectronMomentumWeightedSignal",
                     "Conversion Electron Converted Electron MomentumWeightedSignal", 15, elow, ehigh);
  _cELConvertedElectronCosThetaWeighted =
    tfs->make<TH1D>( "cELConvertedElectronCosThetaWeighted",
                     "Conversion Electron Converted Electron CosThetaWeighted", 200, -1., 1.);
  _cELConvertedElectronCosThetaWeightedSignal =
    tfs->make<TH1D>( "cELConvertedElectronCosThetaWeightedSignal",
                     "Conversion Electron Converted Electron CosThetaWeightedSignal", 200, -1., 1.);


  //
  // and now demand a tracker hit -- but not weighted, since tracker hit is like an acceptance is like a weight
  _cELConvertedElectronMomentumHitTracker =
    tfs->make<TH1D>( "cELConvertedElectronMomentumHitTracker",
                     "Conversion Electron Converted Electron Momentum and Hit Tracker", 200, 0., 200.);
  //
  // this plot is the original momentum of electrons that end up in signal region after energy loss so ehigh
  // needs to be bigger
  _cELConvertedElectronMomentumSignalHitTracker =
    tfs->make<TH1D>( "cELConvertedElectronMomentumSignalHitTracker",
                     "Conversion Electron Converted Electron Momentum Signal and Hit Tracker", 15, elow, 110.);
  _cELConvertedElectronCosThetaHitTracker =
    tfs->make<TH1D>( "cELConvertedElectronCosThetaHitTracker",
                     "Conversion Electron Converted Electron CosTheta Hit Tracker", 200, -1., 1.);
  _cELConvertedElectronCosThetaSignalHitTracker =
    tfs->make<TH1D>( "cELConvertedElectronCosThetaSignalHitTracker",
                     "Conversion Electron Converted Electron CosTheta Signal Hit Tracker", 200, -1., 1.);
  _cELConvertedElectronCosThetaAtTracker=
    tfs->make<TH1D>( "cELConvertedElectronCosThetaAtTracker",
                     "Conversion Electron Converted Electron CosTheta At Tracker", 200, -1., 1.);
  _cELConvertedElectronCosThetaSignalAtTracker =
    tfs->make<TH1D>( "_cELConvertedElectronCosThetaSignalAtTracker",
                     "Conversion Electron Converted Electron CosTheta Signal At Tracker", 200, -1., 1.);
  _cELConvertedElectronMomentumAtTracker =
    tfs->make<TH1D>("cELConvertedElectronMomentumAtTracker","Electron Momentum At Tracker", 200,0.,200.);
  _cELConvertedElectronMomentumSignalAtTracker =
    tfs->make<TH1D>("cELConvertedElectronMomentumSignalAtTracker","Electron Momentum Signal At Tracker", 15,elow,ehigh);
  _cELMomentumLostBeforeTracker =
    tfs->make<TH1D>("_cELMomentumLostBeforeTracker","Electron Momentum Lost Before Tracker", 100,0.,10.);



  _numberOfHitStraws =
    tfs->make<TH1D>("numberOfHitStraws"," Number Of Hit Straws", 100,0.,100.);

  _cELZDiff =
    tfs->make<TH1D>("cELZDiff","Length of Track", 100,0.,5000.);


  //
  // set up dE/dx spectrum

} //bookEventHistos


const double CEL::energyLossSpectrum(const double e){
  //
  // all this from fitting energy loss plotted elsewhere in this job and fitting
  double loss = 0.;
  if (e <= 10.)
    { loss = 4691.*TMath::Landau(e,1.29229,.38467);}// fudge numbers from histo in denom
  if (e >= 2.8 && e <= 20) {loss = TMath::Exp(6.28+ -0.264*(e)); }
  if (e<0){
    throw cet::exception("RANGE")
      << "Nonsense energy in CEL_plugin.cc="
      << e
      << "\n";
  }


  // normalize loss histo to unity (above depended on statistics of particular job)
  return loss/225968.;
}

// Compute a binned representation of the energy loss spectrum.
std::vector<double> CEL::binnedEnergyLossSpectrum(){

  // Sanity check.
  if (_dEdXnbins <= 0) {
    throw cet::exception("RANGE")
      << "Nonsense CEL_plugin.nbins requested="
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
} // CEL::binnedEnergyLossSpectrum


DEFINE_ART_MODULE(CEL);
