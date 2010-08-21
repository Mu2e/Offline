//
// An EDProducer Module that checks radiative pi decays
//
// $Id: RadiativePi_plugin.cc,v 1.1 2010/08/21 17:29:54 rhbob Exp $
// $Author: rhbob $ 
// $Date: 2010/08/21 17:29:54 $
//
// Original author R. Bernstein
//

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>

// Framework includes.
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Services/interface/TFileService.h"
#include "FWCore/Framework/interface/TFileDirectory.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

// Root includes.
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TF1.h"
#include "TSpectrum.h"
#include "TSpectrum2.h"
#include "TSpectrum3.h"
#include "TMultiLayerPerceptron.h"
#include "TMLPAnalyzer.h"


// Mu2e includes.
#include "ToyDP/inc/SimParticleCollection.hh"
#include "GeneralUtilities/inc/RootNameTitleHelper.hh"
#include "GeneralUtilities/inc/pow.hh"
#include "ToyDP/inc/PhysicalVolumeInfoCollection.hh"
#include "Mu2eUtilities/inc/PDGCode.hh"

//CLHEP includes
#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/Randomize.h"

using namespace std;
using CLHEP::Hep3Vector;
using CLHEP::RandPoisson;
using namespace mu2e;
namespace mu2e {

  //--------------------------------------------------------------------
  //
  // 

  class RadiativePi : public edm::EDProducer{
  public:
    RadiativePi(edm::ParameterSet const& pset):
      _maxFullPrint(pset.getUntrackedParameter<int>("maxFullPrint",10)),
      _nAnalyzed(0),
      _messageCategory("ToyHitInfo"){}

    virtual ~RadiativePi() {}

    virtual void beginJob(edm::EventSetup const&);
    virtual void endJob();

    virtual void beginRun(edm::Run const &r, 
                          edm::EventSetup const& eSetup );

    virtual void beginLuminosityBlock(edm::LuminosityBlock const& lblock, 
                                      edm::EventSetup const&);
 
    // This is called for each event.
    void produce(edm::Event& e, edm::EventSetup const&);


  private:

    TH1D* _piCaptureConvertedElectronMomentum;
    TH1D* _piCaptureConvertedElectronMomentumSignal;
    TH1D* _piCaptureConvertedPositronMomentum;
    TH1D* _piCaptureConvertedPositronMomentumSignal;

    TH1D* _piCaptureConvertedElectronCosTheta;
    TH1D* _piCaptureConvertedElectronCosThetaSignal;
    TH1D* _piCaptureConvertedPositronCosTheta;
    TH1D* _piCaptureConvertedPositronCosThetaSignal;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    // Number of events analyzed.
    int _nAnalyzed;


    // A category for the error logger.
    const std::string _messageCategory;

    //name of the module that created the hits to be used
    const std::string _hitCreatorName;


    void RadiativePi::bookEventHistos(double const elow, double const ehigh);
    void RadiativePi::fillEventHistos();
  };


  void RadiativePi::beginJob(edm::EventSetup const& ){

    // Get access to the TFile service.
    //    edm::Service<edm::TFileService> tfs;

  }

  void RadiativePi::endJob(){
  }



  void RadiativePi::beginRun(edm::Run const& run,
			     edm::EventSetup const& eSetup ){
  }

  void RadiativePi::beginLuminosityBlock(edm::LuminosityBlock const& lblock,
					 edm::EventSetup const&){
  }


  void RadiativePi::produce(edm::Event& event, edm::EventSetup const&) {

    static int ncalls(0);
    ++ncalls;

    // Maintain a counter for number of events seen.
    ++_nAnalyzed;

    //    cout << "ncalls = " << ncalls << endl; assert(2==1);


    // Book histogram on the first call regardless
    double elow = 100.;
    double ehigh = 110.;

    if ( ncalls == 1) bookEventHistos(elow,ehigh);

    edm::LogVerbatim log(_messageCategory);
    log << "RadiativePi event #: " 
        << event.id();

    // 
    // start looking through SimParticles
    edm::Handle<SimParticleCollection> simParticles;
    event.getByType(simParticles);

    // Handle to information about G4 physical volumes.
    edm::Handle<PhysicalVolumeInfoCollection> volumes;
    event.getRun().getByType(volumes);

    // Some files might not have the SimParticle and volume information.
    bool haveSimPart = ( simParticles.isValid() && volumes.isValid() );

    // Other files might have empty collections.
    if ( haveSimPart ){
      haveSimPart = !(simParticles->empty() || volumes->empty());
    }

    if (haveSimPart){
      for (uint32_t ithPart = 0; ithPart < simParticles->size(); ++ithPart){
        SimParticle const& sim = simParticles->at(ithPart);
	PhysicalVolumeInfo const& startVol = volumes->at(sim.startVolumeIndex());
	//	cout << " volumename = " << startVol.name() << endl;
	//
	// check three things:  (1) the mother is the original photon, (2) you're an e+ or e-, and (3) the photon converts in the foil
	if ( sim.parentId() == 0 && startVol.name() == "TargetFoil_" ){
	//	if ( sim.parentId() == 0 && startVol.name() == "ToyDSCoil" ){
	  if (sim.pdgId() == PDGCode::e_minus) {
	  CLHEP::HepLorentzVector electronMomentum = sim.startMomentum();
	  //	  cout << "just a test: " << electronMomentum.invariantMass() << endl;
	  double momentum = sqrt(pow(electronMomentum.e(),2) - pow(electronMomentum.invariantMass(),2));
	  _piCaptureConvertedElectronMomentum->Fill(momentum);
	  _piCaptureConvertedElectronMomentumSignal->Fill(momentum);
	  _piCaptureConvertedElectronCosTheta->Fill( electronMomentum.cosTheta() );
	  if (momentum > elow && momentum < ehigh) _piCaptureConvertedElectronCosThetaSignal->Fill(electronMomentum.cosTheta() );
	  
	}
	  if (sim.pdgId() == PDGCode::e_plus) {
	  CLHEP::HepLorentzVector electronMomentum = sim.startMomentum();
	  //	  cout << "just a test: " << electronMomentum.invariantMass() << endl;
	  double momentum = sqrt(pow(electronMomentum.e(),2) - pow(electronMomentum.invariantMass(),2));
	  _piCaptureConvertedPositronMomentum->Fill(momentum);
	  _piCaptureConvertedPositronMomentumSignal->Fill(momentum);
	  _piCaptureConvertedPositronCosTheta->Fill(electronMomentum.cosTheta() );
	  if (momentum > elow && momentum < ehigh) _piCaptureConvertedPositronCosThetaSignal->Fill(electronMomentum.cosTheta() );
	  }
	}
      }
    }
  }
  void RadiativePi::bookEventHistos(double const elow, double const ehigh)
  {        
    //    cout << "booking histos" << endl; assert(2==1);
    edm::Service<edm::TFileService> tfs;
    _piCaptureConvertedElectronMomentum = 
      tfs->make<TH1D>( "piCaptureConvertedElectronMomentum",
		       "Pi Capture Converted Electron Momentum", 200, 0., 200.);
    _piCaptureConvertedElectronMomentumSignal = 
      tfs->make<TH1D>( "piCaptureConvertedElectronMomentumSignal",
		       "Pi Capture Converted Electron MomentumSignal", 20, elow, ehigh);
    _piCaptureConvertedPositronMomentum = 
      tfs->make<TH1D>( "piCaptureConvertedPositronMomentum",
		       "Pi Capture Converted Positron Momentum", 200, 0., 200.);
    _piCaptureConvertedPositronMomentumSignal = 
      tfs->make<TH1D>( "piCaptureConvertedPositronMomentumSignal",
		       "Pi Capture Converted Positron MomentumSignal", 20, elow, ehigh);

    _piCaptureConvertedElectronCosTheta = 
      tfs->make<TH1D>( "piCaptureConvertedElectronCosTheta",
		       "Pi Capture Converted Electron CosTheta", 200, -1., 1.);
    _piCaptureConvertedElectronCosThetaSignal = 
      tfs->make<TH1D>( "piCaptureConvertedElectronCosThetaSignal",
		       "Pi Capture Converted Electron CosThetaSignal", 200, -1., 1.);
    _piCaptureConvertedPositronCosTheta = 
      tfs->make<TH1D>( "piCaptureConvertedPositronCosTheta",
		       "Pi Capture Converted Positron CosTheta", 200, -1., 1.);
    _piCaptureConvertedPositronCosThetaSignal = 
      tfs->make<TH1D>( "piCaptureConvertedPositronCosThetaSignal",
		       "Pi Capture Converted Positron CosThetaSignal", 200, -1., 1.);



  } //bookEventHistos


} // namespace mu2e


DEFINE_FWK_MODULE(RadiativePi);
