//
//
// Simulate the complete process of the nuclear capture of muons by aluminum atoms
// which results in protons, neutrons and photons
//
//
// $Id: NuclearCaptureGun.cc,v 1.24 2013/05/31 20:04:27 gandr Exp $
// $Author: gandr $
// $Date: 2013/05/31 20:04:27 $
//
// Original author Gianni Onorato
//
//

// C++ includes.
#include <iostream>

// Framework includes
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "ConditionsService/inc/DAQParams.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "EventGenerator/inc/NuclearCaptureGun.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "GeneralUtilities/inc/safeSqrt.hh"
#include "StoppingTargetGeom/inc/StoppingTarget.hh"
#include "StoppingTargetGeom/inc/zBinningForFoils.hh"

// Other external includes.
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Units/PhysicalConstants.h"

//ROOT Includes
#include "TH1D.h"
#include "TMath.h"

using namespace std;

static const double spectrumEndPointNeutron = 100.;

namespace mu2e {

  NuclearCaptureGun::NuclearCaptureGun( art::Run& run, const SimpleConfig& config ):

    // Base class.
    GeneratorBase(),

    // Configurable parameters
    _mean(config.getDouble("nuclearCaptureGun.mean",1.)),
    _protonMean(config.getDouble("nuclearCaptureGun.protonMean",0.1)),
    _neutronMean(config.getDouble("nuclearCaptureGun.neutronMean",1.5)),
    _photonMean(config.getDouble("nuclearCaptureGun.photonMean",2.)),
    _protonElow(config.getDouble("nuclearCaptureGun.protonElow",0.)),
    _protonEhi(config.getDouble("nuclearCaptureGun.protonEhi",.300)),
    _neutronElow(config.getDouble("nuclearCaptureGun.neutronElow",0.)),
    _neutronEhi(config.getDouble("nuclearCaptureGun.neutronEhi",spectrumEndPointNeutron)),
    _photonElow(config.getDouble("nuclearCaptureGun.photonElow",0.)),
    _photonEhi(config.getDouble("nuclearCaptureGun.photonEhi",.100)),
    _czmin(config.getDouble("nuclearCaptureGun.czmin",  -1.)),
    _czmax(config.getDouble("nuclearCaptureGun.czmax",  1.)),
    _phimin(config.getDouble("nuclearCaptureGun.phimin", 0. )),
    _phimax(config.getDouble("nuclearCaptureGun.phimax", CLHEP::twopi )),
    _PStoDSDelay(config.getBool("nuclearCaptureGun.PStoDSDelay", true)),
    _pPulseDelay(config.getBool("nuclearCaptureGun.pPulseDelay", false)),
    _pPulseShift(config.getDouble("nuclearCaptureGun.pPulseShift", 0)),
    _timeFolding(config.getBool("FoilParticleGenerator.foldingTimeOption", true)),
    _nProtonBins(config.getInt("nuclearCaptureGun.nProtonBins",1000)),
    _nNeutronBins(evaluateNeutronBins()),
    _nPhotonBins(config.getInt("nuclearCaptureGun.nPhotonBins",1000)),
    _STfname(config.getString("FoilParticleGenerator.STfilename","ExampleDataFiles/StoppedMuons/stoppedMuons_02.txt")),
    _nToSkip (config.getInt("nuclearCaptureGun.nToSkip",0)),
    _doHistograms(config.getBool("nuclearCaptureGun.doHistograms",true)),
    // Initialize random number distributions; getEngine comes from the base class.
    _randPoissonQ( getEngine(), std::abs(_mean) ),
    _randPoissonP( getEngine(), std::abs(_protonMean) ),
    _randPoissonN( getEngine(), std::abs(_neutronMean) ),
    _randPoissonG( getEngine(), std::abs(_photonMean) ),
    _randomUnitSphere ( getEngine(), _czmin, _czmax, _phimin, _phimax ),
    _shapeP ( getEngine() , &(binnedEnergySpectrumProton()[0]), _nProtonBins ),
    _shapeN ( getEngine() , &(binnedEnergySpectrumNeutron()[0]), _nNeutronBins ),
    _shapeG ( getEngine() , &(binnedEnergySpectrumPhoton()[0]), _nPhotonBins ),
    _detSys(),
    // Histogram pointers
    _hNuclearCaptureMultiplicity(),
    _hProtonMultiplicity(),
    _hNeutronMultiplicity(),
    _hPhotonMultiplicity(),
    _hProtonKE(),
    _hNeutronKE(),
    _hPhotonKE(),
    _hProtonKEZoom(),
    _hNeutronKEZoom(),
    _hPhotonKEZoom(),
    _hProtonMomentumMeV(),
    _hNeutronMomentumMeV(),
    _hPhotonMomentumMeV(),
    _hProtonzPosition(),
    _hNeutronzPosition(),
    _hPhotonzPosition(),
    _hProtonCz(),
    _hNeutronCz(),
    _hPhotonCz(),
    _hProtonPhi(),
    _hNeutronPhi(),
    _hPhotonPhi(),
    _hProtonTime(),
    _hNeutronTime(),
    _hPhotonTime(),
    _hProtonMudelay(),
    _hProtonPulsedelay(),
    _hNeutronMudelay(),
    _hNeutronPulsedelay(),
    _hPhotonMudelay(),
    _hPhotonPulsedelay(){

    // About the ConditionsService:
    // The argument to the constructor is ignored for now.  It will be a
    // data base key.  There is a second argument that I have let take its
    // default value of "current"; it will be used to specify a version number.
    ConditionsHandle<AcceleratorParams> accPar("ignored");
    ConditionsHandle<DAQParams>         daqPar("ignored");
    GlobalConstantsHandle<ParticleDataTable> pdt;

    //Set particle mass
    const HepPDT::ParticleData& p_data     = pdt->particle(PDGCode::p_plus).ref();
    const HepPDT::ParticleData& n_data     = pdt->particle(PDGCode::n0).ref();
    _pMass = p_data.mass().value();
    _nMass = n_data.mass().value();

    // Default values for the start and end of the live window.
    // Can be overriden by the run-time config; see below.
    _tmin = daqPar->t0;
    _tmax = accPar->deBuncherPeriod;

    _tmin = config.getDouble("nuclearCaptureGun.tmin",  _tmin );
    _tmax = config.getDouble("nuclearCaptureGun.tmax",  _tmax );

    _detSys = &*GeomHandle<DetectorSystem>();

    // Book histograms.
    if ( _doHistograms ){
      GeomHandle<StoppingTarget> target;
      Binning bins = zBinningForFoils(*target,7);

      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory tfdir  = tfs->mkdir( "NuclearCaptureGun" );

      _hNuclearCaptureMultiplicity = tfdir.make<TH1D>( "hNucCaptMultiplicity", "Multiplicity of Nuclear Capture Events",   20,     0,     20  );
      _hProtonMultiplicity = tfdir.make<TH1D>( "hProtonMultiplicity", "Proton Multiplicity",                20,     0,     20  );
      _hProtonKE           = tfdir.make<TH1D>( "hProtonKE",           "Proton Kinetic Energy",              50, _protonElow,   _protonEhi  );
      _hProtonMomentumMeV  = tfdir.make<TH1D>( "hProtonMomentumMeV",  "Proton Momentum in MeV",             50, _protonElow,   _protonEhi  );
      _hProtonKEZoom       = tfdir.make<TH1D>( "hProtonEZoom",        "Proton Kinetic Energy (zoom)",      200, _protonElow,   _protonEhi  );
      _hProtonzPosition    = tfdir.make<TH1D>( "hProtonzPosition",    "Proton z Position (Tracker Coord)", bins.nbins(), bins.low(), bins.high());
      _hProtonCz           = tfdir.make<TH1D>( "hProtoncz",           "Proton cos(theta)",                 100,    -1.,     1. );
      _hProtonPhi          = tfdir.make<TH1D>( "hProtonphi",          "Proton azimuth",                    100,  -M_PI,  M_PI  );
      _hProtonTime         = tfdir.make<TH1D>( "hProtontime",         "Proton time ",                      210,   -200.,  3000. );
      _hProtonMudelay      = tfdir.make<TH1D>( "hProtonmudelay",      "Proton delay due to muons arriving at ST;(ns)", 300, 0., 2000. );
      _hProtonPulsedelay   = tfdir.make<TH1D>( "hProtonpulsedelay",   "Proton delay due to the proton pulse;(ns)", 60, 0., 300. );

      _hNeutronMultiplicity = tfdir.make<TH1D>( "hNeutronMultiplicity", "Neutron Multiplicity",                20,     0,     20  );
      _hNeutronKE           = tfdir.make<TH1D>( "hNeutronKE",           "Neutron Kinetic Energy",              50, _neutronElow,   _neutronEhi  );
      _hNeutronMomentumMeV  = tfdir.make<TH1D>( "hNeutronMomentumMeV",  "Neutron Momentum in MeV",             50, _neutronElow,   _neutronEhi  );
      _hNeutronKEZoom       = tfdir.make<TH1D>( "hNeutronEZoom",        "Neutron Kinetic Energy (zoom)",      200, _neutronElow,   _neutronEhi  );
      _hNeutronzPosition    = tfdir.make<TH1D>( "hNeutronzPosition",    "Neutron z Position (Tracker Coord)", bins.nbins(), bins.low(), bins.high());
      _hNeutronCz           = tfdir.make<TH1D>( "hNeutroncz",           "Neutron cos(theta)",                 100,    -1.,     1. );
      _hNeutronPhi          = tfdir.make<TH1D>( "hNeutronphi",          "Neutron azimuth",                    100,  -M_PI,  M_PI  );
      _hNeutronTime         = tfdir.make<TH1D>( "hNeutrontime",         "Neutron time ",                      210,   -200.,  3000. );
      _hNeutronMudelay      = tfdir.make<TH1D>( "hNeutronmudelay",      "Neutron delay due to muons arriving at ST;(ns)", 300, 0., 2000. );
      _hNeutronPulsedelay   = tfdir.make<TH1D>( "hNeutronpulsedelay",   "Neutron delay due to the proton pulse;(ns)", 60, 0., 300. );

      _hPhotonMultiplicity = tfdir.make<TH1D>( "hPhotonMultiplicity", "Photon Multiplicity",                20,     0,     20  );
      _hPhotonKE           = tfdir.make<TH1D>( "hPhotonKE",           "Photon Kinetic Energy",              50, _photonElow,   _photonEhi  );
      _hPhotonMomentumMeV  = tfdir.make<TH1D>( "hPhotonMomentumMeV",  "Photon Momentum in MeV",             50, _photonElow,   _photonEhi  );
      _hPhotonKEZoom       = tfdir.make<TH1D>( "hPhotonEZoom",        "Photon Kinetic Energy (zoom)",      200, _photonElow,   _photonEhi  );
      _hPhotonzPosition    = tfdir.make<TH1D>( "hPhotonzPosition",    "Photon z Position (Tracker Coord)", bins.nbins(), bins.low(), bins.high());
      _hPhotonCz           = tfdir.make<TH1D>( "hPhotoncz",           "Photon cos(theta)",                 100,    -1.,     1. );
      _hPhotonPhi          = tfdir.make<TH1D>( "hPhotonphi",          "Photon azimuth",                    100,  -M_PI,  M_PI  );
      _hPhotonTime         = tfdir.make<TH1D>( "hPhotontime",         "Photon time ",                      210,   -200.,  3000. );
      _hPhotonMudelay      = tfdir.make<TH1D>( "hPhotonmudelay",      "Photon delay due to muons arriving at ST;(ns)", 300, 0., 2000. );
      _hPhotonPulsedelay   = tfdir.make<TH1D>( "hPhotonpulsedelay",   "Photon delay due to the proton pulse;(ns)", 60, 0., 300. );
    }

    _fGenerator = unique_ptr<FoilParticleGenerator>(new FoilParticleGenerator( getEngine(), _tmin, _tmax,
                                                                             FoilParticleGenerator::muonFileInputFoil,
                                                                             FoilParticleGenerator::muonFileInputPos,
                                                                             FoilParticleGenerator::negExp,
                                                                             _PStoDSDelay,
                                                                             _pPulseDelay,
                                                                             _pPulseShift,
                                                                             _STfname,
                                                                             _nToSkip));
  }

  NuclearCaptureGun::~NuclearCaptureGun(){
  }




  void NuclearCaptureGun::generate( GenParticleCollection& genParts ){

    // Choose the number of nuclear capture event to generate this event.
    long n = (_mean < 0 ? static_cast<long>(-_mean): _randPoissonQ.fire());

    if ( _doHistograms ) {
        _hNuclearCaptureMultiplicity->Fill(n);
    }

    //Loop over nuclear captures to generate

    long totalP(0), totalN(0), totalG(0);

    for (int i=0; i<n; ++i) {

      //Pick up position and momentum
      CLHEP::Hep3Vector pos(0,0,0);
      double time;
      _fGenerator->generatePositionAndTime(pos, time, _timeFolding);

      const CLHEP::Hep3Vector detPos(_detSys->toDetector(pos));

      //Define the number of protons, neutrons and photons to generate

      long np = (_protonMean < 0 ? static_cast<long>(-_protonMean) : _randPoissonP.fire());
      long nn = (_neutronMean < 0 ? static_cast<long>(-_neutronMean) : _randPoissonN.fire());
      long ng = (_photonMean < 0 ? static_cast<long>(-_photonMean) : _randPoissonG.fire());

      totalP += np;
      totalN += nn;
      totalG += ng;


      //Loop over particles to generate

      for (int ip = 0; ip<np; ++ip) {

          //Pick up energy
        double eKine = _protonElow + _shapeP.fire() * ( _protonEhi - _protonElow );
        double e   = eKine + _pMass;

        //Pick up momentum vector

        _p = safeSqrt(e*e - _pMass*_pMass);
        CLHEP::Hep3Vector p3 = _randomUnitSphere.fire(_p);

        //Set Four-momentum
        CLHEP::HepLorentzVector mom(p3,e);

        // Add the particle to  the list.
        genParts.push_back( GenParticle(PDGCode::p_plus, GenId::nuclearCaptureGun, pos, mom, time));

        // Fill histograms.
        if ( _doHistograms) {
          _hProtonKE->Fill( eKine );
          _hProtonKEZoom->Fill( eKine );
          _hProtonMomentumMeV->Fill( _p );
          _hProtonzPosition->Fill( detPos.z() );
          _hProtonCz->Fill( mom.cosTheta() );
          _hProtonPhi->Fill( mom.phi() );
          _hProtonTime->Fill( time );
          _hProtonMudelay   ->Fill(_fGenerator->muDelay());
          _hProtonPulsedelay->Fill(_fGenerator->pulseDelay());
        }
      } // end of loop on protons

      for (int ineu = 0; ineu<nn; ++ineu) {

        //Pick up energy
        double eKine = _neutronElow + _shapeN.fire() * ( _neutronEhi - _neutronElow );
        double e   = eKine + _nMass;

        //Pick up momentum vector

        _p = safeSqrt(e*e - _nMass*_nMass);
        CLHEP::Hep3Vector p3 = _randomUnitSphere.fire(_p);

        //Set Four-momentum
        CLHEP::HepLorentzVector mom(p3,e);

        // Add the particle to  the list.
        genParts.push_back( GenParticle(PDGCode::n0, GenId::nuclearCaptureGun, pos, mom, time));

        // Fill histograms.
        if ( _doHistograms) {
          _hNeutronKE->Fill( eKine );
          _hNeutronKEZoom->Fill( eKine );
          _hNeutronMomentumMeV->Fill( _p );
          _hNeutronzPosition->Fill( detPos.z() );
          _hNeutronCz->Fill( mom.cosTheta() );
          _hNeutronPhi->Fill( mom.phi() );
          _hNeutronTime->Fill( time );
          _hNeutronMudelay   ->Fill(_fGenerator->muDelay());
          _hNeutronPulsedelay->Fill(_fGenerator->pulseDelay());
        }
      } // end of loop on neutrons

      for (int ig = 0; ig<ng; ++ig) {

        //Pick up energy
        double e   = _photonElow + _shapeG.fire() * ( _photonEhi - _photonElow );

        //Pick up momentum vector

        CLHEP::Hep3Vector p3 = _randomUnitSphere.fire(e);

        //Set Four-momentum
        CLHEP::HepLorentzVector mom(p3,e);

        // Add the particle to  the list.
        genParts.push_back( GenParticle(PDGCode::gamma, GenId::nuclearCaptureGun, pos, mom, time));

        // Fill histograms.
        if ( _doHistograms) {
          _hPhotonKE->Fill( e );
          _hPhotonKEZoom->Fill( e );
          _hPhotonMomentumMeV->Fill( e );
          _hPhotonzPosition->Fill( detPos.z() );
          _hPhotonCz->Fill( mom.cosTheta() );
          _hPhotonPhi->Fill( mom.phi() );
          _hPhotonTime->Fill( time );
          _hPhotonMudelay   ->Fill(_fGenerator->muDelay());
          _hPhotonPulsedelay->Fill(_fGenerator->pulseDelay());
        }
      } // end of loop on photons
    }


    if ( _doHistograms) {
      _hProtonMultiplicity->Fill(totalP);
      _hNeutronMultiplicity->Fill(totalN);
      _hPhotonMultiplicity->Fill(totalG);
    }

  } // end generate

  // Evaluate the number of bins for the neutron spectrum.
  // It has to be compatible with the spectrum table (0.5 KeV)

  int NuclearCaptureGun::evaluateNeutronBins() {

    return (int) ( (_neutronEhi - _neutronElow) / 0.0005 );

  }


  // Energy spectrum of the electron from DIO.
  // Input energy in MeV
  double NuclearCaptureGun::protonEnergySpectrum( double e )
  {

    //taken from GMC
    //
    //   Ed Hungerford  Houston University May 17 1999
    //   Rashid Djilkibaev New York University (modified) May 18 1999
    //
    //   e - proton kinetic energy (MeV)
    //   p - proton Momentum (MeV/c)
    //
    //   Generates a proton spectrum similar to that observed in
    //   u capture in Si.  JEPT 33(1971)11 and PRL 20(1967)569

    //these numbers are in MeV!!!!
    static const double emn = 1.4; // replacing par1 from GMC
    static const double par2 = 1.3279;
    static const double par3=17844.0;
    static const double par4=.32218;
    static const double par5=100.;
    static const double par6=10.014;
    static const double par7=1050.;
    static const double par8=5.103;

    double spectrumWeight;

    if (e >= 20)
      {
        spectrumWeight=par5*TMath::Exp(-(e-20.)/par6);
      }

    else if(e >= 8.0 && e <= 20.0)
      {
        spectrumWeight=par7*exp(-(e-8.)/par8);
      }
    else if (e > emn)
      {
        double xw=(1.-emn/e);
        double xu=TMath::Power(xw,par2);
        double xv=par3*TMath::Exp(-par4*e);
        spectrumWeight=xv*xu;
      }
    else
      {
        spectrumWeight = 0.;
      }
    return spectrumWeight;
  }

  // Compute a binned representation of the energy spectrum of the proton.
  std::vector<double> NuclearCaptureGun::binnedEnergySpectrumProton(){

    // Sanity check.
    if (_nProtonBins <= 0) {
      throw cet::exception("RANGE")
        << "Nonsense nbins requested in "
        << "nuclearCaptureGun (proton) = "
        << _nProtonBins
        << "\n";
    }

    // Bin width.
    double dE = (_protonEhi - _protonElow) / _nProtonBins;

    // Vector to hold the binned representation of the energy spectrum.
    std::vector<double> spectrum;
    spectrum.reserve(_nProtonBins);

    for (int ib=0; ib<_nProtonBins; ib++) {
      double x = _protonElow+(ib+0.5) * dE;
      spectrum.push_back(protonEnergySpectrum(x));
    }

    return spectrum;
  } // NuclearCaptureGun::binnedEnergySpectrumProton


  // Compute a binned representation of the energy spectrum of the neutron.
  //Energy in MeV
  std::vector<double> NuclearCaptureGun::binnedEnergySpectrumNeutron(){

    //The neutron spectum is taken from MARS simulation

    vector<double> neutronSpectrum;
    mu2e::ConfigFileLookupPolicy spectrumFileName;
    string NeutronFileFIP = spectrumFileName("ConditionsService/data/neutronSpectrum.txt");
    fstream infile(NeutronFileFIP.c_str(), ios::in);
    if (infile.is_open()) {
      double en, val;
      int i=0;
      while (i!=_nNeutronBins) {
        if (infile.eof()) break;
        infile >> en >> val;
        if (en >= _neutronElow && en <= _neutronEhi) {
          neutronSpectrum.push_back(val);
          i++;
        }
      }
    } else {
      cout << "No file associated for the ejected neutron spectrum" << endl;
      for (int i=0; i < _nNeutronBins; ++i) {
        neutronSpectrum.push_back(1);
      }
    }


    return neutronSpectrum;
  } // end binnedEnergySpectrumNeutron

  // Compute a binned representation of the energy spectrum of the photon.
  // Energy in MeV
  std::vector<double> NuclearCaptureGun::binnedEnergySpectrumPhoton(){

    vector<double> photonSpectrum;
    for (int i=0; i < _nNeutronBins; ++i) {
      photonSpectrum.push_back(1);
    }

    return photonSpectrum;
  } // end binnedEnergySpectrumPhoton


} // namespace mu2e
