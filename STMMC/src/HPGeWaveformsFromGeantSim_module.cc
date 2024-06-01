// Simulate the electronics response of the HPGe detector. Simulates the pulse height, decay tail, and ADC digitization. Generate one STMWaveformDigi per micropulse.
// Based heavily on example provided in docDb43617 (C. Alvarez-Garcia)
// Pawel Plesniak, 2024
// TODO - some steppoints are dropped as they don't appear to be inside the HPGe volume - figure out why the position is outside the volume.

// stdlib includes
#include <string>
#include <iostream>
#include <bits/stdc++.h>
#include <random>
#include <utility>
#include <algorithm>

// art includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/EDProducer.h"
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

namespace mu2e {
  //================================================================
  class HPGeWaveformsFromGeantSim : public art::EDProducer {
  public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    struct Config {
      // art variable
      fhicl::Atom<std::string> StepPointMCsTag{ Name("StepPointMCsTag"), Comment("InputTag for StepPointMCs")};

      //ADC variables
      fhicl::Atom<uint32_t> fADC{ Name("fADC"), Comment("ADC operating frequency [MHz}")};
      fhicl::OptionalAtom<double> tStep{ Name("tStep"), Comment("Time step for STMWaveformDigis [ns]")};
      fhicl::Atom<double> ADCToEnergy {Name("EnergyPerADCBin"), Comment("ADC energy calibration [keV/bin]")};
      fhicl::OptionalAtom<uint16_t> ADCOffset {Name("ADCOffset"), Comment("Offset from start of ADC charge collection curve. Only use for demo/plotting [ADC counts]")};
      fhicl::Atom<double> noiseSD {Name("NoiseSD"), Comment("Standard deviation of ADC noise [mV]. Set this to 0.0 for the ideal case.")};

      // Physics variables
      fhicl::Atom<double> risingEdgeDecayConstant{ Name("risingEdgeDecayConstant"), Comment("Rising edge decay time [us]")};

      // IO variable
      fhicl::OptionalAtom<bool> verbose {Name("verbose"), Comment("Flag for printing output")};
    };
    using Parameters = art::EDProducer::Table<Config>;
    explicit HPGeWaveformsFromGeantSim(const Parameters& config);

  private:
    void produce(art::Event& event) override;
    void endJob() override;
    void depositCharge(const StepPointMC& step);
    void decayCharge();
    void addNoise();
    void adcs();

    // fhicl variables must be declared first
    std::string _stepPointMCsTag; // Tag of StepPointMCs in STMDet
    uint32_t _fADC; // ADC sampling frequency [MHz]
    double _tStep; // Time step used for simulating the ADC values [ns]
    double _ADCToEnergy; // Calibration of bin width to energy [keV/bin]
    uint16_t _ADCOffset; // Number of ADCs to insert before pulse
    double _noiseSD; // Standard deviation of ADC noise [mV]
    double _risingEdgeDecayConstant; // [us]
    bool _verbose;

    // Define experiment specific constants
    const double _feedbackCapacitance = 1e-12; // [Farads]
    const double _epsilonGe = 2.96; // Energy required to generate an eh pair in Ge at 77K [eV]
    const double _micropulseTime = 1695.0; // [ns]

    // Define physics constants
    const double _e = 1.6e-19; // Electric charge constant [Coulombs]
    const double _driftVelocity = 0.1; // Apprixmate charged particle drift velocity [mm/ns]

    // ADC variables
    double _chargeToADC; // Conversion factor from charge built in capacitor to ADC determined voltage, multiply by this value to get from charge built to ADC voltage.
    uint _nADCs; // Number of ADC values in an event

    // Define Ge crystal properties [mm]
    const double crystalL = 78.5; // Crystal length
    const double crystalR = 36.05; // Crystal radius
    const double crystalHoleL = 64.7; // Crystal hole length not including the hemisphere
    const double crystalHoleR = 5.25; // Crystal hole radius
    const double crystalHoleZStart = crystalL - crystalHoleL; // Starting z position of the crystal hole not including the hemisphere
    const double crystalDirectionGradientCutoff = -crystalHoleZStart/crystalR; // Defines a cone under which points travel to the endcap and not the curved cylinder surface

    // Modelling variables
    double _tMax = 0; // Time of last charge collection [ns]
    double _qMax = 0; // Total charge in capacitor at time of last charge collection: Q(_tMax) = _qMax [Coulombs]
    uint32_t _eventTime = 0; // Time stamp to add to STMWaveformDigi
    std::vector<double> _charge; // Buffer to store charge collected from STMDet StepPointMCs
    std::vector<double> _chargeCollected; // Buffer to store charge collected from STMDet StepPointMCs in the given time step
    std::vector<double> _tmp; // Temporary buffer that will store _chargeCollected over the course of the next event
    std::vector<double> _chargeDecayed; // Buffer to store charge collected that decays over time
    std::vector<int16_t> _adcs; // Buffer for storing the ADC values to put into the STMWaveformDigi
    double hitR = 0; // Hit radial distance [mm]
    double hitZ = 0; // Hit axial distance [mm]
    double R0 = 0; // Hit radial position [mm]
    const double R1 = crystalHoleR; // Distance travelled by electrons [mm]
    double R2 = 0; // Distance travelled by holes [mm]
    double electronTravelDistance = 0, holeTravelDistance = 0; // Drift distances [mm]
    double electronTravelTime = 0, holeTravelTime = 0; // Drift times [ns]
    uint electronTravelTimeSteps = 0, holeTravelTimeSteps = 0; // Drift times [steps]
    uint N_ehPairs = 0; // Number of electron hole pairs
    double trigFactor = 0; // Dimensionless constant used for caluclating distance
    double _decayExp = 0; // Amount of decay with each tStep
    double lastEventEndDecayedCharge = 0;
    CLHEP::Hep3Vector hitPosition; // hit position
    uint32_t eventTime = 0;

    // Offline utilities
    mu2e::STMChannel::enum_type _HPGeChannel = static_cast<mu2e::STMChannel::enum_type>(1);
    STMChannel* _channel = new STMChannel(_HPGeChannel);

    // TODO - want to initialize hpgeEndcapCenterPosition and holeHemisphereCenter as consts here, but errors thrown
    uint64_t stepPointsOutsideHPGe = 0, stepPointsCounter = 0, steps3904 = 0;
    float _ETot = 0.0, _EIon = 0.0;
  };


  //================================================================
  HPGeWaveformsFromGeantSim::HPGeWaveformsFromGeantSim(const Parameters& config )
    : art::EDProducer{config},
    _stepPointMCsTag(config().StepPointMCsTag()),
    _fADC(config().fADC()),
    _ADCToEnergy(config().ADCToEnergy()),
    _noiseSD(config().noiseSD()),
    _risingEdgeDecayConstant(config().risingEdgeDecayConstant())
  {
    produces<STMWaveformDigiCollection>();
    // Assign values to OptionalAtoms
    // _tStep
    if (config().tStep.hasValue()){
      _tStep = *std::move(config().tStep());
      if (_tStep == 0){std::cout << "Time step must be non-zero. Assigning value calculated with fADC." << std::endl;};
    };
    if (!(config().tStep.hasValue()) || (_tStep == 0)){_tStep = 1e3/_fADC;};// 1e3 converts [us] to [ns]

    // ADCOffset
    if (config().ADCOffset.hasValue()){_ADCOffset = *std::move(config().ADCOffset());}
    else{_ADCOffset = 0;};

    // _verbose
    if (config().verbose.hasValue()){_verbose = *std::move(config().verbose());}
    else{_verbose = false;}

    // Define physics parameters to use with the model
    // 1e3 converts keV to eV, multiply by this value to get from charge in capacitor to ADC output voltage
    // Chosen to work in units of fundamental charge _e to avoid requiring larger variables to store parameters
    _chargeToADC = _epsilonGe / (_ADCToEnergy * 1e3);

    // Define the decay amount with each step. 1e3 converts [us] to [ns]
    _decayExp = exp(-_tStep/(_risingEdgeDecayConstant*1e3));

    // Convert noise SD from [mV] to [C]
    _noiseSD = _noiseSD * 1e-3 * _feedbackCapacitance / _e;

    // Determine the number of ADC values in each STMWaveformDigi. Increase the number by one due to truncation. At 320MHz, this will be 543 ADC values per microbunch
    float nADCs = (1695/_tStep) + 1;
    _nADCs = (int) nADCs;
    _adcs.resize(_nADCs, 0);
    _chargeDecayed.resize(_nADCs, 0.);
    _tmp.resize(_nADCs, 0.);

    // Double the size to store all of the charge collections. Expect this to only require a small number of points, factor of 2 chosen for simplicity
    _charge.resize(_nADCs*2, 0.);
    _chargeCollected.resize(_nADCs*2, 0.);
  };

  //================================================================
  void HPGeWaveformsFromGeantSim::endJob()
  {
    std::cout << "StepPointMCsTag " << _stepPointMCsTag << std::endl;
    std::cout << "fADC " << _fADC << std::endl;
    std::cout << "tStep " << _tStep << std::endl;
    std::cout << "ADCToEnergy " << _ADCToEnergy << std::endl;
    std::cout << "ADCOffset " << _ADCOffset << std::endl;
    std::cout << "NoiseSD " << _noiseSD << std::endl;
    std::cout << "Verbose " << _verbose << std::endl;
    std::cout << "N stepPointsOutside HPGe " << stepPointsOutsideHPGe << std::endl;
    std::cout << "N stepPoints " << stepPointsCounter << std::endl;
    std::cout << "_decayExp " << _decayExp << std::endl;
    std::cout << "channel " << &_channel << std::endl;
    std::cout << "channelID " << _channel->id() << "." << std::endl;
    STMChannel::printAll(std::cout);
    std::cout << "stepPointsOutsideHPGe " <<  stepPointsOutsideHPGe << std::endl;
    std::cout << "steps3904 " << steps3904 << std::endl;
    std::cout << "_ETot " << _ETot << std::endl;
    std::cout << "_EIon " << _EIon << std::endl;
  };

  //================================================================
  void HPGeWaveformsFromGeantSim::produce(art::Event& event) {
    // Get the hits in the detector
    art::Handle<StepPointMCCollection> hits_STMDet;
    event.getByLabel(_stepPointMCsTag, hits_STMDet);
    if (!(hits_STMDet.isValid())){std::cout << "Invalid steps. Exiting." << std::endl; return;}

    // Remove me
    std::cout << "Found n steps " << hits_STMDet->size() << std::endl;

    // Clear previous buffer vectors
    std::fill(_chargeDecayed.begin(), _chargeDecayed.end(), 0);
    std::fill(_adcs.begin(), _adcs.end(), 0);

    // Add a collection of charge depositions to _charge
    for(auto& hit : *hits_STMDet){
      _ETot += hit.totalEDep();
      _EIon += hit.ionizingEdep();
      std::cout << "tot " << hit.totalEDep() << std::endl;
      std::cout << "ionizing " << hit.ionizingEdep() << std::endl;
      std::cout << "nonIonizing " << hit.nonIonizingEDep()  << std::endl;
      std::cout << "visible " << hit.visibleEDep() << "\n" << std::endl;
      if (hit.ionizingEdep() != 0)
      {
      depositCharge(hit);
      };
    };

    // Decay all of the collected charges
    decayCharge();

    // Add noise
    addNoise();

    // Update the parameter containing the last charge value from the last event
    lastEventEndDecayedCharge = _chargeDecayed[_nADCs];

    // Convert the charge deposition to ADC voltage output.
    for (uint i = 0; i < _nADCs; i++)
    {
      _adcs[i] = (int16_t) std::round(_chargeDecayed[i]*_chargeToADC);
      std::cout << i << " " << _charge[i] << " " << _chargeDecayed[i] << " " << _adcs[i] << std::endl;
    };

    // Update the parameters to carry over to the next event
    _tmp.clear();
    _tmp.assign(_chargeCollected.begin()+_nADCs, _chargeCollected.end());
    _chargeCollected.clear();
    _chargeCollected.assign(_tmp.begin(), _tmp.end());
    _chargeCollected.insert(_chargeCollected.end(), _nADCs, 0.);

    // Define the time from POT to EWM. Remove the 200ns from eventTime definition.
    /*
    art::Handle<ProtonBunchTimeMC> pbtmcHandle;
    event.getByLabel(_pbtmcTag, pbtmcHandle);
    const ProtonBunchTimeMC& pbtmc(*pbtmcHandle);
    _pbtimemc = pbtmc.pbtime_;
    */

    // Create the STMWaveformDigi and insert all the relevant attributes
    eventTime = ( event.id().event() * _micropulseTime + 200) / 25; // 25ns is 40MHz system clock period.
    STMWaveformDigi _waveformDigi(eventTime, _adcs);
    std::unique_ptr<STMWaveformDigiCollection> outputDigis(new STMWaveformDigiCollection);
    outputDigis->emplace_back(_waveformDigi);

    // Add the STMWaveformDigi to the event
    event.put(std::move(outputDigis));

    std::cout << "\n\n\n\n\n\n Total E dep is " << _ETot << std::endl;
    _ETot=0;

    return;
  };

  //==============================================
  void HPGeWaveformsFromGeantSim::depositCharge(const StepPointMC& step)
  {
    // Define variables that couldn't be constructed in the class constructor
    const CLHEP::Hep3Vector hpgeEndcapCenterPosition(-3946.06,0.0,40671.3); // Crystal position in Mu2e co-ordinate system, defined as the correct quantity
    const CLHEP::Hep3Vector holeHemisphereCenter(0.0, 0.0, crystalHoleZStart); // Crystal hole position in local crystal co-ordinates

    // Calcuate the number of eh pairs from the ionizing energy deposition.
    N_ehPairs = step.ionizingEdep()*1e6/_epsilonGe; // 1e6 converts MeV to eV.

    hitPosition = step.position();
    // Only take the StepPoinMCs in the HPGe detector. STMDet is both sensitive volume of both the HPGe and LaBr.
    if (hitPosition.x() > -3904){
      steps3904++;
      return;
    }

    // Tranform the co-ordinate system to be a cylinder in the +z direction
    if (_verbose) {std::cout << "Position " << hitPosition << std::endl;};
    // Shift the position to a local cylindrical co-ordinate system in with the center of the endcap at the origin
    hitPosition -= hpgeEndcapCenterPosition;
    if (_verbose) {std::cout << "Position in local reference " << hitPosition << std::endl;};
    // Rotate the cylinder for the axis to point in the +z direction
    hitPosition.rotateY(45.0*CLHEP::degree);
    if (_verbose) {std::cout << "Position after rotation " << hitPosition << '\n' << std::endl;};
    // Redefine the hit direction.
    hitR = hitPosition.perp();
    hitZ = hitPosition.z();

    // Run checks
    if ((hitZ < 0) || (hitZ > crystalL) || (hitR < 0)|| (hitR > crystalR) ) //(hitPosition.x() < 0) || (hitPosition.y() < 0) ||
    {
      std::cout << hitPosition << std::endl;
      stepPointsOutsideHPGe++;
      return;
    }

    // Both electrons and holes will travel radially in all cases, but the model volume topology is different depending on the step point position
    if(hitZ > crystalHoleZStart)
    {
      // Volume is a cylinder
      // Electrons travel to the center
      electronTravelDistance = hitR - crystalHoleR;

      // Holes travel to the curved surface
      holeTravelDistance = crystalR - hitR;

      // Define variables for charge buildup time
      R0 = hitR;
      R2 = crystalR;
    }
    else
    {
      // Volume is a sphere
      // Shift the origin to be on the crystal hole center
      hitPosition -= holeHemisphereCenter;
      if(_verbose) {std::cout << "Position after local transformation " << hitPosition << std::endl;};

      // Electrons travel to the surface of the hole hemisphere
      electronTravelDistance = hitPosition.mag() - crystalHoleR;
      if(_verbose) {std::cout << "Magnitude " << hitPosition.mag() << std::endl;};
      if(_verbose) {std::cout << "holeTravelDistance " << holeTravelDistance << std::endl;};

      // Redefine the variables
      hitR = hitPosition.perp();
      hitZ = hitPosition.z();
      R0 = hitPosition.mag();
      if(_verbose) {std::cout << "hitR " << hitR << std::endl;};
      if(_verbose) {std::cout << "hitZ " << hitZ << std::endl;};
      if(_verbose) {std::cout << "R0 "   << R0   << std::endl;};

      // Hole motion is dependent on the position - either travelling to the curved crystal surface or the front of the crystal
      if (hitZ > (crystalDirectionGradientCutoff*hitR)){
      // Holes travel to the curved surface
      trigFactor = hitR / R0; // trigFactor has the positive solution by construction
      holeTravelDistance = (crystalR / trigFactor) - electronTravelDistance - crystalHoleR;
      }
      else
      {
      // Holes travel to the front
      trigFactor = -hitZ / R0; // trigFactor has the positive solution by construction
      holeTravelDistance = (crystalHoleZStart / trigFactor) - electronTravelDistance - crystalHoleR;
      }
      R2 = electronTravelDistance + holeTravelDistance + crystalHoleR;
    }

    if ((holeTravelDistance < 0)  || (electronTravelDistance < 0))
    {
      std::cout << "holeTravelDistance is " << holeTravelDistance << std::endl;
      std::cout << ", electronTravelDistance is " << electronTravelDistance << std::endl;
      std::cout << ". Both should be positive but one is not. Exiting." << std::endl;
      exit(-1);
    }

    // Calculate the drift times
    electronTravelTime = electronTravelDistance / _driftVelocity;
    holeTravelTime = holeTravelDistance / _driftVelocity;

    // Define parameters required for charge deposition. Constants A and B are defined here for efficiency
    uint tIndex = step.time()/_tStep;
    const uint tIndexStart = tIndex;
    const double A = N_ehPairs/log(R1/R2); // Should be R2/R1, but multipling this number by -1 caused computational issues. This gets the same correct number.
    const double B = _driftVelocity/R0;

    // Scale and shift the charge particle travel time
    electronTravelTimeSteps = tIndex + electronTravelTime/_tStep;
    holeTravelTimeSteps = tIndex + holeTravelTime/_tStep;

    if(_verbose){
      std::cout << "Summary" << std::endl;
      std::cout << "electronTravelTime " << electronTravelTime << std::endl;
      std::cout << "holeTravelTime " << holeTravelTime << std::endl;
      std::cout << "Number of eh pairs " << N_ehPairs << std::endl;
      std::cout << "step time" << step.time() << std::endl;
      std::cout << "step time in index " << step.time()/_tStep << std::endl;
      std::cout << "tIndex " << tIndex << std::endl;
      std::cout << "_tStep " << _tStep << std::endl;
      std::cout << "electronTravelTimeSteps " << electronTravelTimeSteps << std::endl;
      std::cout << "holeTravelTimeSteps " << holeTravelTimeSteps << std::endl;
    };

    // Calculate charge collection when both particles are moving
    while ((tIndex < electronTravelTimeSteps) & (tIndex < holeTravelTimeSteps))
    {
      _charge[tIndex] = A*log((1+B*(tIndex-tIndexStart)*_tStep)/(1-B*(tIndex-tIndexStart)*_tStep));
      tIndex++;
    };

    // Calculate charge build up when only the holes are moving
    if ((tIndex >= electronTravelTimeSteps) & (tIndex < holeTravelTimeSteps))
    {
      while(tIndex <= holeTravelTimeSteps)
      {
      _charge[tIndex] = A * log( (1+B*(tIndex-tIndexStart)*_tStep ) / ( R1/R0 ) );
      tIndex++;
      };
    }
    // Calculate charge build up when only the electrons are moving
    else if ((tIndex < electronTravelTimeSteps) & (tIndex >= holeTravelTimeSteps))
    {
      while(tIndex <= electronTravelTimeSteps)
      {
      _charge[tIndex] = A * log( ( R2/R0 ) / (1-B*(tIndex-tIndexStart)*_tStep) );
      tIndex++;
      };
    }
    else
    {
      std::cout << "Logic error with generating the correct _charge vector. Exiting." << std::endl;
      exit(-1);
    };

    if(tIndexStart == 0)
    {
      _chargeCollected[tIndexStart] += _charge[tIndexStart];
    }

    for (uint i = tIndexStart; i < tIndex; i++)
    {
      _chargeCollected[i] += _charge[i] - _charge[i-1];
    };

    std::fill(_charge.begin(), _charge.end(), 0);
    return;
  };

  //==============================================
  void HPGeWaveformsFromGeantSim::decayCharge()
  {
    _chargeDecayed[0] = lastEventEndDecayedCharge * _decayExp + _chargeCollected[0];
    for (uint t = 1; t < _nADCs; t++)
    {
      _chargeDecayed[t] = _chargeDecayed[t-1] * _decayExp + _chargeCollected[t];
    };
    return;
  };

  // ==============================================
  void HPGeWaveformsFromGeantSim::addNoise()
  {
    if (_noiseSD == 0.)
    {
      return;
    }
    std::default_random_engine _randomGen;
    std::normal_distribution<double> _noiseDistribution(0.0, _noiseSD);
    for (size_t _i = 0; _i < _nADCs; _i++)
    {
      _chargeDecayed[_i] += _noiseDistribution(_randomGen);
    };
    return;
  };
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::HPGeWaveformsFromGeantSim)
