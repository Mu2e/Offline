// Simulates the electronics response of the HPGe detector. Simulates the pulse height, decay tail, and ADC digitization. Generates one STMWaveformDigi per micropulse.
// Model based heavily on example provided in docDb 43617
// See docDb 51487 for full documentation
// Original author: Pawel Plesniak

// stdlib includes
#include <algorithm>
#include <bits/stdc++.h>
#include <cmath>
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
// #include "Offline/DataProducts/inc/STMChannel.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
// #include "Offline/Mu2eUtilities/inc/STMUtils.hh"
// #include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/STMWaveformDigi.hh"
// #include "Offline/STMConditions/inc/STMEnergyCalib.hh"

// ROOT includes
#include "art_root_io/TFileService.h"
#include "TTree.h"


namespace mu2e {
  class HPGeWaveformsFromStepPointMCs : public art::EDProducer {
  public:
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;
    struct Config {
      fhicl::Atom<art::InputTag> StepPointMCsTagEle{ Name("StepPointMCsTagEle"), Comment("InputTag for StepPointMCs derived from EleBeamCat")};
      fhicl::Atom<art::InputTag> StepPointMCsTagMu{ Name("StepPointMCsTagMu"), Comment("InputTag for StepPointMCs derived from MuBeamCat")};
      fhicl::Atom<art::InputTag> StepPointMCsTag1809{ Name("StepPointMCsTag1809"), Comment("InputTag for StepPointMCs derived from TargetStopsCat")};
      fhicl::Atom<double> fADC{ Name("fADC"), Comment("ADC operating frequency [MHz}")};
      fhicl::Atom<double> ADCToEnergy {Name("EnergyPerADCBin"), Comment("ADC energy calibration [keV/bin]")};
      fhicl::Atom<double> noiseSD {Name("NoiseSD"), Comment("Standard deviation of ADC noise [mV]. Set this to 0.0 for the ideal case.")};
      fhicl::Atom<double> risingEdgeDecayConstant{ Name("risingEdgeDecayConstant"), Comment("Rising edge decay time [us]")};
      fhicl::OptionalAtom<int> microspillBufferLengthCount{ Name("microspillBufferLengthCount"), Comment("Number of microspills to buffer ahead for, in number of microspills")};
      fhicl::OptionalAtom<bool> makeTTree{ Name("makeTTree"), Comment("Controls whether to make the TTree with branches chargeCollected, chargeDecayed, ADC, eventId, time")};
      fhicl::OptionalAtom<double> timeOffset{ Name("timeOffset"), Comment("For debugging, adds the named time offset in [ns], used for testing analysis algorithms")};
      fhicl::OptionalAtom<uint> resetEventNumber{ Name("resetEventNumber"), Comment("Simulates the off-spill period by resetting the inter-event last decayed charge to zero")};
      fhicl::OptionalAtom<uint> verbosityLevel{ Name("verbosityLevel"), Comment("Controls verbosity")};
    };
    using Parameters = art::EDProducer::Table<Config>;
    explicit HPGeWaveformsFromStepPointMCs(const Parameters& conf);
  private:
    void produce(art::Event& event) override;
    void beginJob();
    void depositCharge(const StepPointMC& step);
    void decayCharge();
    void addNoise();
    void digitize();

    // fhicl variables
    art::ProductToken<StepPointMCCollection> StepPointMCsTokenEle, StepPointMCsTokenMu, StepPointMCsToken1809;  // Token of StepPointMCs in STMDet
    double fADC = 0;                                                                                            // ADC sampling frequency [MHz]
    double ADCToEnergy = 0;                                                                                     // Calibration of bin width to energy [keV/bin]
    double noiseSD = 0;                                                                                         // Standard deviation of ADC noise [mV]
    double risingEdgeDecayConstant = 0;                                                                         // [us]
    bool makeTTree = false;                                                                                     // Controls whether an analysis TTree is made
    double timeOffset = 0.0;                                                                                    // Used for debugging [ns]
    uint resetEventNumber = 0;                                                                                  // Event ID at which to reset the charge amplitude
    int verbosityLevel = 0;                                                                                     // How much output to generate

    // Define experiment specific constants
    const double feedbackCapacitance = 1e-12;                                                                   // [Farads]
    const double epsilonGe = 2.96;                                                                              // Energy required to generate an eh pair in Ge at 77K [eV]
    const double micropulseTime = 1695.0;                                                                       // [ns]

    // Define physics constants
    const double _e = 1.602176634e-19;                                                                          // Electric charge constant [Coulombs]
    const double electronDriftVelocity = 0.08;                                                                  // Apprixmate charged particle drift velocity [mm/ns]
    const double holeDriftVelocity = 0.06;                                                                      // Apprixmate charged particle drift velocity [mm/ns]

    // ADC variables
    double chargeToADC = 0;                                                                                     // Conversion factor from charge built in capacitor to ADC determined voltage, multiply by this value to get from charge built to ADC voltage.
    uint nADCs = 0;                                                                                             // Number of ADC values in an event
    const int16_t ADCMax = static_cast<int16_t>((-1 * std::pow(2, 15)) + 1);                                    // Maximum ADC value, power is 15 not 16 as using int16_t not uint16_t
    double ADC = 0;                                                                                             // iterator variable
    uint32_t eventTimeBuffer = 0;                                                                               // Multiple of event ids to store

    // Define Ge crystal properties [mm]
    const double crystalCentreX = -3973.81;                                                                     // Crystal centre x position [mm]
    const double crystalCentreY = 0;                                                                            // Crystal centre y position [mm]
    const double crystalCentreZ = 40699.1;                                                                      // Crystal centre z position [mm]
    CLHEP::Hep3Vector crystalCentrePosition;                                                                    // Crystal centre position vector
    // TODO - want to initialize hpgeEndcapCenterPosition and holeHemisphereCenter as consts here, but errors thrown
    CLHEP::Hep3Vector hitPosition;
    const double crystalL = 78.5;                                                                               // Crystal length [mm]
    const double crystalR = 36.05;                                                                              // Crystal radius [mm]
    const double crystalHoleL = 64.7;                                                                           // Crystal hole length not including the hemisphere [mm]
    const double crystalHoleR = 5.25;                                                                           // Crystal hole radius [mm]
    const double crystalHoleZStart = crystalL - crystalHoleL;                                                   // Starting z position of the crystal hole not including the hemisphere [mm]
    const double crystalDirectionGradientCutoff = -crystalHoleZStart/crystalR;                                  // Defines a cone under which points travel to the endcap and not the curved cylinder surface
    const double stepPositionTolerance = 0.1;                                                                   // Adjusts for resolution of applying rotation
    const double maxR = crystalR + stepPositionTolerance;                                                       // Crystal radius including tolderance
    const double maxZ = crystalL + stepPositionTolerance;                                                       // Crystal z including tolerance

    // Modelling variables
    double hitR = 0;                                                                                            // Hit radial distance [mm]
    double hitZ = 0;                                                                                            // Hit axial distance [mm]
    double R0 = 0;                                                                                              // Hit radial position [mm]
    const double R1 = crystalHoleR;                                                                             // Distance travelled by electrons [mm]
    double R2 = 0;                                                                                              // Distance travelled by holes [mm]

    double trigFactor = 0;                                                                                      // Dimensionless constant used for caluclating distance
    int32_t N_ehPairs = 0;                                                                                      // Number of electron hole pairs
    uint32_t eventTime = 0;                                                                                     // Time stamp to add to STMWaveformDigi [ADC clock ticks]

    double tADC = 0;                                                                                            // Time step used for simulating the ADC values [ns]
    double electronTravelDistance = 0, holeTravelDistance = 0;                                                  // Drift distances [mm]
    double electronTravelTime = 0, holeTravelTime = 0;                                                          // Drift times [ns]
    uint32_t electronTravelTimeSteps = 0, holeTravelTimeSteps = 0;                                              // Drift times [steps]
    double decayExp = 0;                                                                                        // Amount of decay with each tADC
    double lastEventEndDecayedCharge = 0;                                                                       // Carry over for starting new microspill waveforms [charge carrier pairs]
    int microspillBufferLengthCount = 0;                                                                        // Buffer to store the charge deposits that are allocated to this event but happen after the microspill ends e.g. 844keV
    const int defaultMicrospillBufferLengthCount = 2;                                                           // Default value for the microspill buffer length

    // TTree and storage variables
    TTree* ttree;                                                                                               // ttree variable
    double chargeCollected = 0, chargeDecayed = 0;                                                              // used to fill the ttree
    uint eventId = 0;                                                                                           // used to fill the ttree
    uint32_t time = 0;                                                                                          // used to fill the ttree
    int16_t ttreeADC = 0;                                                                                       // used to fill the ttree

    // Data storage vectors
    std::vector<double> _charge;                                                                                // Buffer to store charge collected from STMDet StepPointMCs
    std::vector<double> _chargeCollected;                                                                       // Buffer to store charge collected from STMDet StepPointMCs in the given time step
    std::vector<double> _chargeDecayed;                                                                         // Buffer to store charge collected that decays over time
    std::vector<double> _chargeCarryOver;                                                                       // Temporary buffer that will store _chargeCollected over the course of the next event
    std::vector<int16_t> _adcs;                                                                                 // Buffer for storing the ADC values to put into the STMWaveformDigi

    // Offline utilities
    // TODO: include the prodition to get the sampling frequency
    // mu2e::STMChannel::enum_type _HPGeChannel = static_cast<mu2e::STMChannel::enum_type>(1);
    // STMChannel* _channel = new STMChannel(_HPGeChannel);
    // ProditionsHandle<STMEnergyCalib> _stmEnergyCalib_h;
  };

  HPGeWaveformsFromStepPointMCs::HPGeWaveformsFromStepPointMCs(const Parameters& conf)
    : art::EDProducer{conf},
      StepPointMCsTokenEle(consumes<StepPointMCCollection>(conf().StepPointMCsTagEle())),
      StepPointMCsTokenMu(consumes<StepPointMCCollection>(conf().StepPointMCsTagMu())),
      StepPointMCsToken1809(consumes<StepPointMCCollection>(conf().StepPointMCsTag1809())),
      fADC(conf().fADC()),
      ADCToEnergy(conf().ADCToEnergy()),
      noiseSD(conf().noiseSD()),
      risingEdgeDecayConstant(conf().risingEdgeDecayConstant()) {
        produces<STMWaveformDigiCollection>();
        if (defaultMicrospillBufferLengthCount < 2)
          throw cet::exception("RANGE", "defaultMicrospillBufferLengthCount has to be more than 1\n");

        crystalCentrePosition.set(crystalCentreX, crystalCentreY, crystalCentreZ);
        tADC = 1e3/fADC; // 1e3 converts [us] to [ns] as fADC is in [MHz]

        // Assign optional variables
        microspillBufferLengthCount = conf().microspillBufferLengthCount() ? *(conf().microspillBufferLengthCount()) : defaultMicrospillBufferLengthCount;
        verbosityLevel = conf().verbosityLevel() ? *(conf().verbosityLevel()) : 0;

        // Determine the number of ADC values in each STMWaveformDigi. Increase the number by one due to truncation. At 320MHz, this will be 543 ADC values per microbunch
        double _nADCs = (micropulseTime/tADC) + 1;
        nADCs = (int) _nADCs;
        _charge.insert(_charge.begin(), nADCs * microspillBufferLengthCount, 0.);
        _chargeCollected.insert(_chargeCollected.begin(), nADCs * microspillBufferLengthCount, 0.);
        _chargeCarryOver.insert(_chargeCarryOver.begin(), nADCs * (microspillBufferLengthCount - 1), 0.);
        _chargeDecayed.insert(_chargeDecayed.begin(), nADCs, 0.);
        _adcs.insert(_adcs.begin(), nADCs, 0);

        // Define physics parameters to use with the model
        // Define the decay amount with each step. 1e3 converts [us] to [ns]
        decayExp = exp(-tADC/(risingEdgeDecayConstant*1e3));

        // Convert noise SD from [mV] to [C], division by _e to work in units of charge
        noiseSD *= (1e-3 * feedbackCapacitance/_e);

        // 1e3 converts keV to eV, multiply by this value to get from charge in capacitor to ADC output voltage
        // Chosen to work in units of fundamental charge _e to avoid storing very small double values
        chargeToADC = epsilonGe / (ADCToEnergy * 1e3);

        // Assign the approrpiate time offset
        timeOffset = conf().timeOffset() ? *(conf().timeOffset()) : 0.0;

        // Assign TTrees
        makeTTree = conf().makeTTree() ? *(conf().makeTTree()) : false;
        if (makeTTree) {
          art::ServiceHandle<art::TFileService> tfs;
          ttree = tfs->make<TTree>("ttree", "MakeHPGeWaveformsFromStepPointMCs ttree");
          ttree->Branch("chargeCollected", &chargeCollected, "chargeCollected/D");
          ttree->Branch("chargeDecayed", &chargeDecayed, "chargeDecayed/D");
          ttree->Branch("ADC", &ttreeADC, "ADC/S");
          ttree->Branch("eventId", &eventId, "eventId/i");
          ttree->Branch("time", &time, "time/i");
        };

        resetEventNumber = conf().resetEventNumber() ? *(conf().resetEventNumber()) : 0;
      };

  void HPGeWaveformsFromStepPointMCs::beginJob() {
    if (verbosityLevel) {
      std::cout << "STM HPGe digitization parameters" << std::endl;
      std::cout << "\tInput parameters" << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "fAD [MHz]"                            << fADC                                     << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "EnergyPerADCBin [keV/bin]"            << ADCToEnergy                              << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "NoiseSD [mV]"                         << noiseSD /(1e-3 * feedbackCapacitance/_e) << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "risingEdgeDecayConstant [us]"         << risingEdgeDecayConstant                  << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "microspillBufferLengthCount"          << microspillBufferLengthCount              << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "makeTTree"                            << makeTTree                                << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "timeOffset [ns]"                      << timeOffset                               << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "resetEventNumber"                     << resetEventNumber                         << std::endl;
      std::cout << "\tDerived parameters: " << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "tADC [ns]"                            << tADC                                     << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "nADCs"                                << nADCs                                    << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "NoiseSD [charge carriers]"            << noiseSD                                  << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "chargeToADC [bin/charge carrier]"     << chargeToADC                              << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "Voltage range [V]"                    << "[+1, -1]"                               << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "Voltage range used [V]"               << "[0, -1]"                                << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "Voltage range used [bins]"            << "[0, " << ADCMax << "]"                  << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "Voltage range used [charge carriers]" << "[0, " << ADCMax/chargeToADC << "]"      << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "Voltage range used [C]"               << "[0, " << ADCMax * _e/chargeToADC << "]" << std::endl;
      std::cout << std::left << "\t\t" << std::setw(60) << "Energy range [keV]"                   << "[0, " << -1*ADCMax * ADCToEnergy << "]" << std::endl;
      std::cout << std::endl; // buffer line
    };
  };

  void HPGeWaveformsFromStepPointMCs::produce(art::Event& event) {
    eventId = event.id().event();
    // Get the hits in the detector
    std::vector<StepPointMC> StepsEle = event.getProduct(StepPointMCsTokenEle);
    std::vector<StepPointMC> StepsMu = event.getProduct(StepPointMCsTokenMu);
    std::vector<StepPointMC> Steps1809 = event.getProduct(StepPointMCsToken1809);

    // Add a collection of charge depositions to _charge
    for(const StepPointMC& step : StepsEle){
      if (step.ionizingEdep() != 0)
          depositCharge(step);
    };
    for(const StepPointMC& step : StepsMu){
      if (step.ionizingEdep() != 0)
          depositCharge(step);
    };
    for(const StepPointMC& step : Steps1809){
      if (step.ionizingEdep() != 0)
          depositCharge(step);
    };

    // Decay all of the collected charges
    decayCharge();

    // Update the last event decayed charge before the noise so the noise isn't added twice
    if (resetEventNumber != 0 && eventId == resetEventNumber) {
      lastEventEndDecayedCharge = 0;
      std::fill(_chargeCarryOver.begin(), _chargeCarryOver.end(), 0);
    }
    else
      lastEventEndDecayedCharge = _chargeDecayed.back();

    // Add preamplifier electronics noise with SD defined in noiseSD
    addNoise();

    // Digitize the waveform
    digitize();

    // Validate the digitized waveforms
    for (auto j : _adcs) {
      if (j > 1000)
        throw cet::exception("LogicError", "ADC values too high!");
    };

    // Simulation takes the POT time as t = 0, and has sequential microspills (events). The trigger time offset is not used here, left as a TODO
    // Create the STMWaveformDigi and insert all the relevant attributes
    // TODO - this only keeps accurate time if the sampling is 320MHz. Needs to be rewritten to work for other times
    eventTimeBuffer = eventId % 5;
    if (eventTimeBuffer == 0 || (eventId % 3) == 0)
      eventTime += nADCs + 1;
    else
      eventTime += nADCs;
    STMWaveformDigi _waveformDigi(eventTime, _adcs);
    std::unique_ptr<STMWaveformDigiCollection> outputDigis(new STMWaveformDigiCollection);
    outputDigis->emplace_back(_waveformDigi);

    // Make the ttree if appropriate
    if (makeTTree) {
      time = _waveformDigi.trigTimeOffset() - 1;
      for (uint i = 0; i < nADCs; i++) {
        chargeCollected = _chargeCollected[i];
        chargeDecayed = _chargeDecayed[i];
        ttreeADC = _adcs[i];
        time++;
        ttree->Fill();
      };
    };

    // Update the parameters to carry over to the next event
    _chargeCarryOver.clear();
    _chargeCarryOver.assign(_chargeCollected.begin() + nADCs, _chargeCollected.end());
    _chargeCollected.clear();
    _chargeCollected.assign(_chargeCarryOver.begin(), _chargeCarryOver.end());
    _chargeCollected.insert(_chargeCollected.end(), nADCs, 0);
    // Clear previous buffer vectors
    std::fill(_chargeDecayed.begin(), _chargeDecayed.end(), 0);
    std::fill(_adcs.begin(), _adcs.end(), 0);

    // Add the STMWaveformDigi to the event
    event.put(std::move(outputDigis));
    return;
  };

  void HPGeWaveformsFromStepPointMCs::depositCharge(const StepPointMC& step) {
    // Define variables that couldn't be constructed in the class constructor
    const CLHEP::Hep3Vector holeHemisphereCenter(0.0, 0.0, crystalHoleZStart); // Crystal hole position in local crystal co-ordinates

    hitPosition = step.position();
    // Only take the StepPoinMCs in the HPGe detector. STMDet is both sensitive volume of both the HPGe and LaBr.
    if (hitPosition.x() > -3904)
      return;
    // If the time is outside the buffer time, skip it
    if (step.time() > (microspillBufferLengthCount * micropulseTime))
      return;

    // Tranform the co-ordinate system to be a cylinder in the +z direction
    // Shift the position to a local cylindrical co-ordinate system with the center of the crystal at the origin
    hitPosition -= crystalCentrePosition;
    // Rotate the cylinder for the axis to point in the +z direction
    hitPosition.rotateY(45.0*CLHEP::degree);
    // Shift the crystal so the front of the crystal is the start
    hitPosition.setZ(hitPosition.z() + (crystalL/2));

    // Redefine the hit direction.
    hitR = hitPosition.perp();
    hitZ = hitPosition.z();

    // Run checks
    if (hitZ < -stepPositionTolerance)
      throw cet::exception("LogicError") << "Step not inside HPGe detector, z position negative: " << hitZ << "\n";
    if (hitZ > maxZ)
      throw cet::exception("LogicError") << "Step not inside HPGe detector, z position greater than the crystal length: " << hitZ << ", should be in range [0, " << maxZ << "]\n";
    if (hitR > maxR)
      throw cet::exception("LogicError") << "Step not inside HPGe detector, radius outside of crystal: " << hitR << ", should be in range [0, " << maxR << "].\n";

    // Both electrons and holes will travel radially in all cases, but the model volume topology is different depending on the step point position
    if (hitZ > crystalHoleZStart) { // Volume is a cylinder
      electronTravelDistance = hitR + stepPositionTolerance - crystalHoleR; // Electrons travel to the center
      // Holes travel to the curved surface
      holeTravelDistance = crystalR + stepPositionTolerance - hitR;
      // Define variables for charge buildup time
      R0 = hitR;
      R2 = crystalR + stepPositionTolerance * 2;
    }
    else { // Volume is a sphere
      // Shift the origin to be on the crystal hole center
      hitPosition -= holeHemisphereCenter;
      // Redefine the variables
      hitR = hitPosition.perp();
      hitZ = hitPosition.z();
      R0 = hitPosition.mag();
      // Electrons travel to the surface of the hole hemisphere
      electronTravelDistance = R0 + stepPositionTolerance - crystalHoleR;

      // Hole motion is dependent on the position - either travelling to the curved crystal surface or the front of the crystal
      // Note - hole motion is modelled to always be radial. Required as field here is effectively radial.
      if (hitZ > (crystalDirectionGradientCutoff * hitR)) { // Holes travel to the curved surface
        trigFactor = hitR / R0; // trigFactor is positive solution by construction
        holeTravelDistance = (crystalR / trigFactor) - R0 + stepPositionTolerance;
      }
      else { // Holes travel to the endcap
        trigFactor = -hitZ / R0; // trigFactor is positive solution by construction
        holeTravelDistance = (crystalHoleZStart / trigFactor) - R0 + stepPositionTolerance;
      };
      R2 = electronTravelDistance + holeTravelDistance + crystalHoleR + stepPositionTolerance * 2;
    };

    if (holeTravelDistance < 0 || electronTravelDistance < 0)
      throw cet::exception("LogicError") << "Electron (" << electronTravelDistance << ") and hole (" << holeTravelDistance << ") travelling distances should both be positve.\nPosition found at " << step.position() << "(" << hitPosition << ")\n";

    // Calculate the drift times
    electronTravelTime = electronTravelDistance / electronDriftVelocity;
    holeTravelTime = holeTravelDistance / holeDriftVelocity;

    // Calcuate the number of eh pairs from the ionizing energy deposition.
    N_ehPairs = -1.0 * step.ionizingEdep() * 1e6 / epsilonGe; // 1e6 converts MeV to eV. -1.0 as this is a decreasing peak

    // Define parameters required for charge deposition. Constants A and B are defined here for code brevity
    uint tIndex = (step.time() + timeOffset) / tADC, tIndexStart = tIndex;
    const double A = N_ehPairs / log(R2/R1);
    const double Be = electronDriftVelocity / R0;
    const double Bh = holeDriftVelocity / R0;

    // Scale and shift the charge particle travel time
    electronTravelTimeSteps = tIndex + electronTravelTime/tADC;
    holeTravelTimeSteps = tIndex + holeTravelTime/tADC;

    // Calculate charge collection when both particles are moving
    while (tIndex <= electronTravelTimeSteps && tIndex <= holeTravelTimeSteps) {
      _charge[tIndex] = A*log((1 + Bh * (tIndex - tIndexStart) * tADC) / (1 - Be * (tIndex - tIndexStart) * tADC));
      tIndex++;
    };

    // Calculate charge build up when only the holes are moving
    if (tIndex >= electronTravelTimeSteps && tIndex < holeTravelTimeSteps) {
      while(tIndex < holeTravelTimeSteps) {
        _charge[tIndex] = A * log((1 + Bh * (tIndex - tIndexStart) * tADC ) / (R1 / R0));
        tIndex++;
      };
    }
    // Calculate charge build up when only the electrons are moving
    else if (tIndex < electronTravelTimeSteps && tIndex >= holeTravelTimeSteps) {
      while(tIndex < electronTravelTimeSteps) {
        _charge[tIndex] = A * log((R2 / R0) / (1 - Be * (tIndex - tIndexStart) * tADC));
        tIndex++;
      };
    };

    // Allocate the rest of the charge for one more entry for the continuity
    _charge[tIndex] = N_ehPairs;
    tIndex++;

    // Update _chargeCollected. First case is treated separately as there is no charge deposited in the previous step
    if(tIndexStart == 0) {
      _chargeCollected[tIndexStart] += _charge[tIndexStart];
      tIndexStart++;
    };
    for (uint i = tIndexStart; i < tIndex; i++)
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
    for (uint i = 0; i < nADCs; i++) {
      ADC = _chargeDecayed[i] * chargeToADC;
      _adcs[i] = ADC > ADCMax ? static_cast<int16_t>(std::round(ADC)) : ADCMax;
    };
    return;
  };
}; // namespace mu2e

DEFINE_ART_MODULE(mu2e::HPGeWaveformsFromStepPointMCs)
