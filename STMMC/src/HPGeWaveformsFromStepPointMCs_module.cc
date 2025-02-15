// Simulates the electronics response of the HPGe detector. Simulates the pulse height, decay tail, and ADC digitization. Generates one STMWaveformDigi per micropulse.
// Model based heavily on example provided in docDb 43617
// See docDb 51487 for full documentation
// Input parameters
//   - StepPointMCsTag - StepPointMCs in the HPGe detector
//   - fADC - ADC sampling frequency, in MHz
//   - ADCToEnergy - ADC energy calibration, in keV/bin
//   - noiseSD - Standard deviation of the electronics noise, in keV
//   - risingEdgeDecayConstant - rising edge decay time constant, in us
//   - HPGeCrystalEndcapCentreX - x of the HPGe crystal endcap centre
//   - HPGeCrystalEndcapCentreY - y of the HPGe crystal endcap centre
//   - HPGeCrystalEndcapCentreZ - z of the HPGe crystal endcap centre
// Original author: Pawel Plesniak

// stdlib includes
#include <algorithm>
#include <bits/stdc++.h>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <utility>

// art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Units/PhysicalConstants.h"

// exception handling
#include "cetlib_except/exception.h"

// fhicl includes
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"

// message handling
#include "messagefacility/MessageLogger/MessageLogger.h"

// Offline includes
#include "Offline/DataProducts/inc/STMChannel.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/Mu2eUtilities/inc/STMUtils.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/STMWaveformDigi.hh"
#include "Offline/STMConditions/inc/STMEnergyCalib.hh"

// ROOT includes
#include "art_root_io/TFileService.h"
#include "TTree.h"


namespace mu2e {
  class HPGeWaveformsFromStepPointMCs : public art::EDProducer {
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
      fhicl::Atom<int> microspillBufferLengthCount{ Name("microspillBufferLengthCount"), Comment("Number of microspills to buffer ahead for")};
    };
    using Parameters = art::EDProducer::Table<Config>;
    explicit HPGeWaveformsFromStepPointMCs(const Parameters& conf);
  private:
    void produce(art::Event& event) override;
    void depositCharge(const StepPointMC& step);
    void decayCharge();
    void addNoise();
    void digitize();

    art::ProductToken<StepPointMCCollection> StepPointMCsToken; // Token of StepPointMCs in STMDet
    uint32_t fADC = 0; // ADC sampling frequency [MHz]
    double tStep = 0; // Time step used for simulating the ADC values [ns]
    double ADCToEnergy = 0; // Calibration of bin width to energy [keV/bin]
    double noiseSD = 0; // Standard deviation of ADC noise [mV]
    double risingEdgeDecayConstant = 0; // [us]
    double HPGeCrystalEndcapCentreX = 0.0, HPGeCrystalEndcapCentreY = 0.0, HPGeCrystalEndcapCentreZ = 0.0; // Endcap centre variables

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
    double masterClockTickPeriod = 25.0; //ns

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

    double trigFactor = 0; // Dimensionless constant used for caluclating distance
    int32_t N_ehPairs = 0; // Number of electron hole pairs
    uint32_t eventTime = 0; // Time stamp to add to STMWaveformDigi

    double electronTravelDistance = 0, holeTravelDistance = 0; // Drift distances [mm]
    double electronTravelTime = 0, holeTravelTime = 0; // Drift times [ns]
    uint32_t electronTravelTimeSteps = 0, holeTravelTimeSteps = 0; // Drift times [steps]
    double decayExp = 0; // Amount of decay with each tStep

    // TODO - want to initialize hpgeEndcapCenterPosition and holeHemisphereCenter as consts here, but errors thrown
    CLHEP::Hep3Vector hitPosition; // hit position
    std::vector<double> _charge; // Buffer to store charge collected from STMDet StepPointMCs
    std::vector<double> _chargeCollected; // Buffer to store charge collected from STMDet StepPointMCs in the given time step
    std::vector<double> _chargeCarryOver; // Temporary buffer that will store _chargeCollected over the course of the next event
    std::vector<double> _chargeDecayed; // Buffer to store charge collected that decays over time
    std::vector<int16_t> _adcs; // Buffer for storing the ADC values to put into the STMWaveformDigi
    int microspillBufferLengthCount;
    const int defaultMicrospillBufferLengthCount = 2;
    double lastEventEndDecayedCharge = 0; // Carry over for starting new microspill waveforms

    // Offline utilities
    mu2e::STMChannel::enum_type _HPGeChannel = static_cast<mu2e::STMChannel::enum_type>(1);
    STMChannel* _channel = new STMChannel(_HPGeChannel);
    ProditionsHandle<STMEnergyCalib> _stmEnergyCalib_h;

  };

  HPGeWaveformsFromStepPointMCs::HPGeWaveformsFromStepPointMCs(const Parameters& conf )
    : art::EDProducer{conf},
      StepPointMCsToken(consumes<StepPointMCCollection>(conf().StepPointMCsTag())),
      fADC(conf().fADC()),
      ADCToEnergy(conf().ADCToEnergy()),
      noiseSD(conf().noiseSD()),
      risingEdgeDecayConstant(conf().risingEdgeDecayConstant()),
      HPGeCrystalEndcapCentreX(conf().HPGeCrystalEndcapCentreX()),
      HPGeCrystalEndcapCentreY(conf().HPGeCrystalEndcapCentreY()),
      HPGeCrystalEndcapCentreZ(conf().HPGeCrystalEndcapCentreZ()) {
        produces<STMWaveformDigiCollection>();
        tStep = 1e3/fADC; // 1e3 converts [us] to [ns] as fADC is in [MHz]

        // Determine the number of ADC values in each STMWaveformDigi. Increase the number by one due to truncation. At 320MHz, this will be 543 ADC values per microbunch
        double _nADCs = (micropulseTime/tStep) + 1;
        nADCs = (int) _nADCs;
        _charge.resize(nADCs * microspillBufferLengthCount, 0.);
        _chargeCollected.resize(nADCs * microspillBufferLengthCount, 0.);
        _chargeCarryOver.resize(nADCs * (microspillBufferLengthCount - 1), 0.);
        _chargeDecayed.resize(nADCs, 0.);
        _adcs.resize(nADCs, 0);

        // Define physics parameters to use with the model
        // Define the decay amount with each step. 1e3 converts [us] to [ns]
        decayExp = exp(-tStep/(risingEdgeDecayConstant*1e3));

        // Convert noise SD from [mV] to [C], division by _e to work in units of charge
        noiseSD *= 1e-3 * feedbackCapacitance/_e;

        // 1e3 converts keV to eV, multiply by this value to get from charge in capacitor to ADC output voltage
        // Chosen to work in units of fundamental charge _e to avoid storing very small double values
        chargeToADC = epsilonGe / (ADCToEnergy * 1e3);

        // Assign microspillBufferLengthCount
        microspillBufferLengthCount = conf().microspillBufferLengthCount() ? conf().microspillBufferLengthCount() : defaultMicrospillBufferLengthCount;
      };

  void HPGeWaveformsFromStepPointMCs::produce(art::Event& event) {
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
    _chargeCarryOver.clear();
    _chargeCarryOver.assign(_chargeCollected.begin() + nADCs, _chargeCollected.end());
    _chargeCollected.clear();
    _chargeCollected.assign(_chargeCarryOver.begin(), _chargeCarryOver.begin() + nADCs);
    _chargeCollected.insert(_chargeCollected.end(), nADCs * (microspillBufferLengthCount - 1), 0.);

    // Define the time from POT to EWM. Remove the 200ns from eventTime definition.
    /*
    art::Handle<ProtonBunchTimeMC> pbtmcHandle;
    event.getByLabel(_pbtmcTag, pbtmcHandle);
    const ProtonBunchTimeMC& pbtmc(*pbtmcHandle);
    _pbtimemc = pbtmc.pbtime_;
    */

    // Create the STMWaveformDigi and insert all the relevant attributes
    eventTime = (event.id().event() * micropulseTime + 200.0) / masterClockTickPeriod;
    STMWaveformDigi _waveformDigi(eventTime, _adcs);
    std::unique_ptr<STMWaveformDigiCollection> outputDigis(new STMWaveformDigiCollection);
    outputDigis->emplace_back(_waveformDigi);

    // Add the STMWaveformDigi to the event
    event.put(std::move(outputDigis));
    return;
  };

  void HPGeWaveformsFromStepPointMCs::depositCharge(const StepPointMC& step) {
    // Define variables that couldn't be constructed in the class constructor
    const CLHEP::Hep3Vector hpgeEndcapCenterPosition(HPGeCrystalEndcapCentreX, HPGeCrystalEndcapCentreY, HPGeCrystalEndcapCentreZ); // Crystal position in Mu2e co-ordinate system, defined as the correct quantity
    const CLHEP::Hep3Vector holeHemisphereCenter(0.0, 0.0, crystalHoleZStart); // Crystal hole position in local crystal co-ordinates

    hitPosition = step.position();
    // Only take the StepPoinMCs in the HPGe detector. STMDet is both sensitive volume of both the HPGe and LaBr.
    if (hitPosition.x() > -3904) {
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
      throw cet::exception("LogicError") << "Step not inside HPGe detector, exiting.\n";

    // Both electrons and holes will travel radially in all cases, but the model volume topology is different depending on the step point position
    if (hitZ > crystalHoleZStart) { // Volume is a cylinder
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
      // Note - hole motion is modelled to always be radial. Required as field here is effectively radial.
      if (hitZ > (crystalDirectionGradientCutoff * hitR)) { // Holes travel to the curved surface
        trigFactor = hitR / R0; // trigFactor is positive solution by construction
        holeTravelDistance = (crystalR / trigFactor) - R0;
      }
      else { // Holes travel to the front
        trigFactor = -hitZ / R0; // trigFactor is positive solution by construction
        holeTravelDistance = (crystalHoleZStart / trigFactor) - R0;
      };
      R2 = electronTravelDistance + holeTravelDistance + crystalHoleR;
    };

    if (holeTravelDistance < 0 || electronTravelDistance < 0)
      throw cet::exception("LogicError") << "Electron (" << electronTravelDistance << ") and hole (" << holeTravelDistance << ") travelling distances should both be positve.\nPosition found at " << step.position() << "(" << hitPosition << ")" ;

    // Calculate the drift times
    electronTravelTime = electronTravelDistance / electronDriftVelocity;
    holeTravelTime = holeTravelDistance / holeDriftVelocity;

    // Calcuate the number of eh pairs from the ionizing energy deposition.
    N_ehPairs = -1.0 * step.ionizingEdep() * 1e6 / epsilonGe; // 1e6 converts MeV to eV. -1.0 as this is a decreasing peak

    // Define parameters required for charge deposition. Constants A and B are defined here for code brevity
    uint tIndex = step.time() / tStep, tIndexStart = tIndex;
    const double A = N_ehPairs / log(R2/R1);
    const double Be = electronDriftVelocity / R0;
    const double Bh = holeDriftVelocity / R0;

    // Scale and shift the charge particle travel time
    electronTravelTimeSteps = tIndex + electronTravelTime/tStep;
    holeTravelTimeSteps = tIndex + holeTravelTime/tStep;

    // Calculate charge collection when both particles are moving
    while (tIndex < electronTravelTimeSteps && tIndex < holeTravelTimeSteps) {
      _charge[tIndex] = A*log((1 + Bh * (tIndex - tIndexStart) * tStep) / (1 - Be * (tIndex - tIndexStart) * tStep));
      tIndex++;
    };

    // Calculate charge build up when only the holes are moving
    if (tIndex >= electronTravelTimeSteps && tIndex < holeTravelTimeSteps) {
      while(tIndex <= holeTravelTimeSteps) {
        _charge[tIndex] = A * log((1 + Bh * (tIndex - tIndexStart) * tStep ) / (R1 / R0));
        tIndex++;
      };
    }
    // Calculate charge build up when only the electrons are moving
    else if (tIndex < electronTravelTimeSteps && tIndex >= holeTravelTimeSteps) {
      while(tIndex <= electronTravelTimeSteps) {
        _charge[tIndex] = A * log((R2 / R0) / (1 - Be * (tIndex - tIndexStart) * tStep));
        tIndex++;
      };
    }
    else
      throw cet::exception("LogicError") << "Error with the charged particle travel time.";

    // Allocate the rest of the charge accumulated over the interaction vertex
    for (uint i = tIndex; i < _charge.size(); i++)
      _charge[i] = N_ehPairs;

    // Update _chargeCollected. First case is treated separately as there is no charge deposited in the previous step
    if(tIndexStart == 0) {
      _chargeCollected[tIndexStart] += _charge[tIndexStart];
      tIndexStart++;
    };
    for (uint i = tIndexStart; i < _charge.size(); i++)
      _chargeCollected[i] += (_charge[i] - _charge[i-1]);

    // Clear the charge vector
    std::fill(_charge.begin(), _charge.end(), 0);
    return;
  };

  void HPGeWaveformsFromStepPointMCs::decayCharge() {
    _chargeDecayed[0] = lastEventEndDecayedCharge * decayExp + _chargeCollected[0];
    for (uint t = 1; t < nADCs; t++)
      _chargeDecayed[t] = _chargeDecayed[t-1] * decayExp + _chargeCollected[t];
    return;
  };

  void HPGeWaveformsFromStepPointMCs::addNoise() {
    // If the noise SD is zero, do nothing
    if (noiseSD < std::numeric_limits<double>::epsilon())
      return;

    // Add the noise in multiples of the fundamental charge
    std::default_random_engine _randomGen;
    std::normal_distribution<double> _noiseDistribution(0.0, noiseSD);
    for (size_t _i = 0; _i < nADCs; _i++)
      _chargeDecayed[_i] += _noiseDistribution(_randomGen);
    return;
  };

  void HPGeWaveformsFromStepPointMCs::digitize() {
    // Convert the charge deposition to ADC voltage output.
    for (uint i = 0; i < nADCs; i++)
      _adcs[i] = (int16_t) std::round(_chargeDecayed[i]*chargeToADC);
    return;
  };
}; // namespace mu2e

DEFINE_ART_MODULE(mu2e::HPGeWaveformsFromStepPointMCs)
