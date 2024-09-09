// Simulate the electronics response of the HPGe detector. Simulates the pulse height, decay tail, and ADC digitization. Generate one STMWaveformDigi per micropulse.
// Based heavily on example provided in docDb43617 (C. Alvarez-Garcia)
// Pawel Plesniak, 2024

// stdlib includes
#include <string>
#include <iostream>
#include <bits/stdc++.h>
#include <random>
#include <utility>
#include <algorithm>

// art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"

// fhicl includes
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "canvas/Utilities/InputTag.h"

// Offline includes
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/RecoDataProducts/inc/STMWaveformDigi.hh"
#include "Offline/Mu2eUtilities/inc/STMUtils.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/STMConditions/inc/STMEnergyCalib.hh"
#include "Offline/DataProducts/inc/STMChannel.hh"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Units/PhysicalConstants.h"

// ROOT includes
#include "TH1.h"
#include "TTree.h"

// Exception handling
#include "cetlib_except/exception.h"

namespace mu2e {
  //================================================================
  class HPGeWaveformsFromGeantSim : public art::EDProducer {
  public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    struct Config {
      fhicl::Atom<art::InputTag> StepPointMCsTag{ Name("StepPointMCsTag"), Comment("InputTag for StepPointMCs")};
      fhicl::Atom<uint32_t> fADC{ Name("fADC"), Comment("ADC operating frequency [MHz}")};
      fhicl::Atom<double> ADCToEnergy {Name("EnergyPerADCBin"), Comment("ADC energy calibration [keV/bin]")};
      fhicl::Atom<double> noiseSD {Name("NoiseSD"), Comment("Standard deviation of ADC noise [mV]. Set this to 0.0 for the ideal case.")};
      fhicl::Atom<double> risingEdgeDecayConstant{ Name("risingEdgeDecayConstant"), Comment("Rising edge decay time [us]")};
      fhicl::Atom<double> HPGeCrystalEndcapCentreX{ Name("HPGeCrystalEndcapCentreX"), Comment("HPGe endcap centre x")};
      fhicl::Atom<double> HPGeCrystalEndcapCentreY{ Name("HPGeCrystalEndcapCentreY"), Comment("HPGe endcap centre y")};
      fhicl::Atom<double> HPGeCrystalEndcapCentreZ{ Name("HPGeCrystalEndcapCentreZ"), Comment("HPGe endcap centre z")};
      fhicl::OptionalAtom<bool> verbose {Name("verbose"), Comment("Flag for printing output")};
    };
    using Parameters = art::EDProducer::Table<Config>;
    explicit HPGeWaveformsFromGeantSim(const Parameters& conf);

  private:
    void produce(art::Event& event) override;
    void depositCharge(const StepPointMC& step);
    void decayCharge();
    void addNoise();
    void digitize();

    // fhicl variables must be declared first
    art::ProductToken<StepPointMCCollection> StepPointMCsToken; // Token of StepPointMCs in STMDet
    uint32_t fADC = 0; // ADC sampling frequency [MHz]
    double tStep = 0; // Time step used for simulating the ADC values [ns]
    double ADCToEnergy = 0; // Calibration of bin width to energy [keV/bin]
    double noiseSD = 0; // Standard deviation of ADC noise [mV]
    double risingEdgeDecayConstant = 0; // [us]
    double HPGeCrystalEndcapCentreX = 0.0, HPGeCrystalEndcapCentreY = 0.0, HPGeCrystalEndcapCentreZ = 0.0; // Endcap centre variables
    bool verbose = false;

    // Define experiment specific constants
    const double feedbackCapacitance = 1e-12; // [Farads]
    const double epsilonGe = 2.96; // Energy required to generate an eh pair in Ge at 77K [eV]
    const double micropulseTime = 1695.0; // [ns]

    // Define physics constants
    const double _e = 1.602176634e-19; // Electric charge constant [Coulombs]
    const double electronDriftVelocity = 0.08; // Apprixmate charged particle drift velocity [mm/ns]
    const double holeDriftVelocity = 0.06; // Apprixmate charged particle drift velocity [mm/ns]

    // ADC variables
    double chargeToADC = 0; // Conversion factor from charge built in capacitor to ADC determined voltage, multiply by this value to get from charge built to ADC voltage.
    uint nADCs = 0; // Number of ADC values in an event

    // Define Ge crystal properties [mm]
    const double crystalL = 78.5; // Crystal length
    const double crystalR = 36.05; // Crystal radius
    const double crystalHoleL = 64.7; // Crystal hole length not including the hemisphere
    const double crystalHoleR = 5.25; // Crystal hole radius
    const double crystalHoleZStart = crystalL - crystalHoleL; // Starting z position of the crystal hole not including the hemisphere
    const double crystalDirectionGradientCutoff = -crystalHoleZStart/crystalR; // Defines a cone under which points travel to the endcap and not the curved cylinder surface

    // Modelling variables
    double hitR = 0; // Hit radial distance [mm]
    double hitZ = 0; // Hit axial distance [mm]
    double R0 = 0; // Hit radial position [mm]
    const double R1 = crystalHoleR; // Distance travelled by electrons [mm]
    double R2 = 0; // Distance travelled by holes [mm]

    double electronTravelDistance = 0, holeTravelDistance = 0; // Drift distances [mm]
    double electronTravelTime = 0, holeTravelTime = 0; // Drift times [ns]
    double trigFactor = 0; // Dimensionless constant used for caluclating distance
    double decayExp = 0; // Amount of decay with each tStep
    double lastEventEndDecayedCharge = 0;

    uint32_t electronTravelTimeSteps = 0, holeTravelTimeSteps = 0; // Drift times [steps]
    int32_t N_ehPairs = 0; // Number of electron hole pairs
    uint32_t eventTime = 0; // Time stamp to add to STMWaveformDigi

    CLHEP::Hep3Vector hitPosition; // hit position
    std::vector<double> _charge; // Buffer to store charge collected from STMDet StepPointMCs
    std::vector<double> _chargeCollected; // Buffer to store charge collected from STMDet StepPointMCs in the given time step
    std::vector<double> _tmp; // Temporary buffer that will store _chargeCollected over the course of the next event
    std::vector<double> _chargeDecayed; // Buffer to store charge collected that decays over time
    std::vector<int16_t> _adcs; // Buffer for storing the ADC values to put into the STMWaveformDigi

    // Offline utilities
    mu2e::STMChannel::enum_type _HPGeChannel = static_cast<mu2e::STMChannel::enum_type>(1);
    STMChannel* _channel = new STMChannel(_HPGeChannel);
    ProditionsHandle<STMEnergyCalib> _stmEnergyCalib_h;

    // TODO - want to initialize hpgeEndcapCenterPosition and holeHemisphereCenter as consts here, but errors thrown
  };


  //================================================================
  HPGeWaveformsFromGeantSim::HPGeWaveformsFromGeantSim(const Parameters& conf )
    : art::EDProducer{conf},
      StepPointMCsToken(consumes<StepPointMCCollection>(conf().StepPointMCsTag())),
      fADC(conf().fADC()),
      ADCToEnergy(conf().ADCToEnergy()),
      noiseSD(conf().noiseSD()),
      risingEdgeDecayConstant(conf().risingEdgeDecayConstant()),
      HPGeCrystalEndcapCentreX(conf().HPGeCrystalEndcapCentreX()),
      HPGeCrystalEndcapCentreY(conf().HPGeCrystalEndcapCentreY()),
      HPGeCrystalEndcapCentreZ(conf().HPGeCrystalEndcapCentreZ())
  {
    produces<STMWaveformDigiCollection>();
    // Assign values to OptionalAtoms
    auto _verbose = conf().verbose();
    if(_verbose)verbose = *_verbose;

    tStep = 1e3/fADC; // 1e3 converts [us] to [ns]

    // Define physics parameters to use with the model
    // 1e3 converts keV to eV, multiply by this value to get from charge in capacitor to ADC output voltage
    // Chosen to work in units of fundamental charge _e to avoid requiring larger variables to store parameters
    chargeToADC = epsilonGe / (ADCToEnergy * 1e3);

    // Define the decay amount with each step. 1e3 converts [us] to [ns]
    decayExp = exp(-tStep/(risingEdgeDecayConstant*1e3)); // DELETEME - WAS 1

    // Convert noise SD from [mV] to [C]
    noiseSD = noiseSD * 1e-3 * feedbackCapacitance/_e;

    // Determine the number of ADC values in each STMWaveformDigi. Increase the number by one due to truncation. At 320MHz, this will be 543 ADC values per microbunch
    float _nADCs = (micropulseTime/tStep) + 1;
    nADCs = (int) _nADCs;
    _adcs.resize(nADCs, 0);
    _chargeDecayed.resize(nADCs, 0.);
    _tmp.resize(nADCs, 0.);

    // Double the size to store all of the charge collections. Expect this to only require a small number of points, factor of 2 chosen for simplicity
    _charge.resize(nADCs*2, 0.);
    _chargeCollected.resize(nADCs*2, 0.);
  };

  //================================================================
  void HPGeWaveformsFromGeantSim::produce(art::Event& event) {
    // Get the hits in the detector
    auto const& hits_STMDet = event.getProduct(StepPointMCsToken);
    // Clear previous buffer vectors
    std::fill(_chargeDecayed.begin(), _chargeDecayed.end(), 0);
    std::fill(_adcs.begin(), _adcs.end(), 0);

    // Add a collection of charge depositions to _charge
    for(auto& hit : hits_STMDet){
      if (hit.ionizingEdep() != 0)
          depositCharge(hit);
    };

    // Decay all of the collected charges
    decayCharge();

    // Add preamplifier electronics noise with SD defined in noiseSD
    addNoise();

    // Digitize the waveform
    digitize();

    // Update the parameters to carry over to the next event
    lastEventEndDecayedCharge = _chargeDecayed[nADCs];
    _tmp.clear();
    _tmp.assign(_chargeCollected.begin()+nADCs, _chargeCollected.end());
    _chargeCollected.clear();
    _chargeCollected.assign(_tmp.begin(), _tmp.end());
    _chargeCollected.insert(_chargeCollected.end(), nADCs, 0.);

    // Define the time from POT to EWM. Remove the 200ns from eventTime definition.
    /*
    art::Handle<ProtonBunchTimeMC> pbtmcHandle;
    event.getByLabel(_pbtmcTag, pbtmcHandle);
    const ProtonBunchTimeMC& pbtmc(*pbtmcHandle);
    _pbtimemc = pbtmc.pbtime_;
    */

    // Create the STMWaveformDigi and insert all the relevant attributes
    eventTime = ( event.id().event() * micropulseTime + 200) / 25; // 25ns is 40MHz system clock period.
    eventTime = eventTime + tStep*20; // DEUBGGING
    STMWaveformDigi _waveformDigi(eventTime, _adcs);
    std::unique_ptr<STMWaveformDigiCollection> outputDigis(new STMWaveformDigiCollection);
    outputDigis->emplace_back(_waveformDigi);

    // Add the STMWaveformDigi to the event
    event.put(std::move(outputDigis));
    return;
  };

  //==============================================
  void HPGeWaveformsFromGeantSim::depositCharge(const StepPointMC& step)
  {
    // Define variables that couldn't be constructed in the class constructor
    const CLHEP::Hep3Vector hpgeEndcapCenterPosition(HPGeCrystalEndcapCentreX, HPGeCrystalEndcapCentreY, HPGeCrystalEndcapCentreZ); // Crystal position in Mu2e co-ordinate system, defined as the correct quantity
    const CLHEP::Hep3Vector holeHemisphereCenter(0.0, 0.0, crystalHoleZStart); // Crystal hole position in local crystal co-ordinates

    hitPosition = step.position();
    // Only take the StepPoinMCs in the HPGe detector. STMDet is both sensitive volume of both the HPGe and LaBr.
    if (hitPosition.x() > -3904){
      return;
    }

    // Tranform the co-ordinate system to be a cylinder in the +z direction
    // Shift the position to a local cylindrical co-ordinate system in with the center of the endcap at the origin
    hitPosition -= hpgeEndcapCenterPosition;
    // Rotate the cylinder for the axis to point in the +z direction
    hitPosition.rotateY(45.0*CLHEP::degree);
    // Redefine the hit direction.
    hitR = hitPosition.perp();
    hitZ = hitPosition.z();

    // Run checks
    if ((hitZ < 0) || (hitZ > crystalL) || (hitR < 0) || (hitR > crystalR))
      throw cet::exception("RANGE") << "Step not inside HPGe detector, exiting.\n";

    // Both electrons and holes will travel radially in all cases, but the model volume topology is different depending on the step point position
    if(hitZ > crystalHoleZStart)
    { // Volume is a cylinder
      electronTravelDistance = hitR - crystalHoleR; // Electrons travel to the center
      // Holes travel to the curved surface
      holeTravelDistance = crystalR - hitR;
      // Define variables for charge buildup time
      R0 = hitR;
      R2 = crystalR;
    }
    else
    { // Volume is a sphere
      // Shift the origin to be on the crystal hole center
      hitPosition -= holeHemisphereCenter;
      // Electrons travel to the surface of the hole hemisphere
      electronTravelDistance = hitPosition.mag() - crystalHoleR;
      // Redefine the variables
      hitR = hitPosition.perp();
      hitZ = hitPosition.z();
      R0 = hitPosition.mag();

      // Hole motion is dependent on the position - either travelling to the curved crystal surface or the front of the crystal
      if (hitZ > (crystalDirectionGradientCutoff*hitR)){
      // Holes travel to the curved surface
        trigFactor = hitR / R0; // trigFactor has the positive solution by construction
        holeTravelDistance = (crystalR / trigFactor) - R0;
      }
      else
      {
      // Holes travel to the front
        trigFactor = -hitZ / R0; // trigFactor has the positive solution by construction
        holeTravelDistance = (crystalHoleZStart / trigFactor) - R0;
      }
      R2 = electronTravelDistance + holeTravelDistance + crystalHoleR;
    }

    if ((holeTravelDistance < 0)  || (electronTravelDistance < 0))
    {
      throw cet::exception("RANGE") << "Electron (" << electronTravelDistance << ") and hole (" << holeTravelDistance << ") travelling distances should both be positve.\nPosition found at " << step.position() << "(" << hitPosition << ")" ;
    }

    // Calculate the drift times
    electronTravelTime = electronTravelDistance / electronDriftVelocity;
    holeTravelTime = holeTravelDistance / holeDriftVelocity;

    // Calcuate the number of eh pairs from the ionizing energy deposition.
    N_ehPairs = -1.0*step.ionizingEdep()*1e6/epsilonGe; // 1e6 converts MeV to eV. -1.0 used as every quantity is later scaled by -1.0.

    // Define parameters required for charge deposition. Constants A and B are defined here for code brevity
    uint tIndex = (step.time()/tStep) + 25; // REMOVETHEPLUS25
    uint tIndexStart = tIndex;

    //const double A = N_ehPairs/log(R1/R2); // Should be R2/R1, but multipling this number by -1 caused computational issues. This gets the same correct number.
    const double A = N_ehPairs/(log(R2/R1)); // Should be R2/R1, but multipling this number by -1 caused computational issues. This gets the same correct number.
    const double Be = electronDriftVelocity/R0;
    const double Bh = holeDriftVelocity/R0;

    // Scale and shift the charge particle travel time
    electronTravelTimeSteps = tIndex + electronTravelTime/tStep;
    holeTravelTimeSteps = tIndex + holeTravelTime/tStep;

    // Calculate charge collection when both particles are moving
    while ((tIndex < electronTravelTimeSteps) && (tIndex < holeTravelTimeSteps))
    {
      _charge[tIndex] = A*log((1+Bh*(tIndex-tIndexStart)*tStep)/(1-Be*(tIndex-tIndexStart)*tStep));
      tIndex++;
    };

    // Calculate charge build up when only the holes are moving
    if ((tIndex >= electronTravelTimeSteps) && (tIndex < holeTravelTimeSteps))
    {
      while(tIndex <= holeTravelTimeSteps)
      {
        _charge[tIndex] = A * log( (1+Bh*(tIndex-tIndexStart)*tStep ) / ( R1/R0 ) );
        tIndex++;
      };
    }
    // Calculate charge build up when only the electrons are moving
    else if ((tIndex < electronTravelTimeSteps) && (tIndex >= holeTravelTimeSteps))
    {
      while(tIndex <= electronTravelTimeSteps)
      {
        _charge[tIndex] = A * log( ( R2/R0 ) / (1-Be*(tIndex-tIndexStart)*tStep) );
        tIndex++;
      };
    }
    else
    {
      throw cet::exception("RANGE") << "Error with the charged particle travel time.";
    };

    // Allocate the rest of the charge accumulated over the interaction vertex
    for (uint i = tIndex; i < (2*nADCs); i++)
      _charge[i] = N_ehPairs;

    for (auto i : _charge)
      std::cout << i << ", ";
    std::cout << std::endl;

    // The charge collected is at the start of this event - DELETEME review this
    if(tIndexStart == 0)
    {
      _chargeCollected[tIndexStart] += _charge[tIndexStart];
      tIndexStart++;
    }

    for (uint i = tIndexStart; i < tIndex; i++)
      _chargeCollected[i] += (_charge[i] - _charge[i-1]);

    // Clear the charge vector
    std::fill(_charge.begin(), _charge.end(), 0);

    return;
  };

  //==============================================
  void HPGeWaveformsFromGeantSim::decayCharge()
  {
    _chargeDecayed[0] = lastEventEndDecayedCharge * decayExp + _chargeCollected[0];
    for (uint t = 1; t < nADCs; t++)
    {
      _chargeDecayed[t] = _chargeDecayed[t-1] * decayExp + _chargeCollected[t];
    };
    return;
  };

  // ==============================================
  void HPGeWaveformsFromGeantSim::addNoise()
  {
    if (noiseSD == 0.)
    {
      return;
    }
    std::default_random_engine _randomGen;
    std::normal_distribution<double> _noiseDistribution(0.0, noiseSD);
    for (size_t _i = 0; _i < nADCs; _i++)
    {
      _chargeDecayed[_i] += _noiseDistribution(_randomGen);
    };
    return;
  };
  // ==============================================
  void HPGeWaveformsFromGeantSim::addNoise()
  {
    // Convert the charge deposition to ADC voltage output.
    for (uint i = 0; i < nADCs; i++)
      _adcs[i] = (int16_t) std::round(_chargeDecayed[i]*chargeToADC);
    return;
  };
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::HPGeWaveformsFromGeantSim)
