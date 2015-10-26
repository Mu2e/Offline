
/*

  A plug_in for running a variety of event generators.

  Original author Rob Kutschke

  Eventually this will support a variety of event generators, controllable
  from the run time configuration.  A given call might invoke one or more
  of these generators.

  1) A full featured single particle gun.
  2) Single conversion track, uniformly from the targets.
  3) (Emax-E)**5 DIO model.
  4) Other DIO models.
  5) protons, neutrons, gammas and nuclear fragments from muon capture.
  6) Mockups of CLHEP::pion capture on nuclei and of CLHEP::pion and muon decay in flight.
  I say mock-ups because I see this starting from an known CLHEP::pion and muon
  flux distributions, not by starting from a CLHEP::pion or a muon entering
  the DS.
  7) Simplified models of cosmics.

  At present I expect that the highest fidelity generation of cosmics will be
  done by running an external generator and then reading "events" from the output
  of that generator.  Perhaps the merge will be done in this module, perhaps
  it will be done in a separate module?

*/

// Mu2e includes.
#include "ConfigTools/inc/SimpleConfig.hh"
#include "ConfigTools/inc/requireUniqueKey.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/G4BeamlineInfoCollection.hh"

// Particular generators that this code knows about.
#include "EventGenerator/inc/ConversionGun.hh"
#include "EventGenerator/inc/StoppedMuonGun.hh"
#include "EventGenerator/inc/CosmicDYB.hh"
#include "EventGenerator/inc/DecayInOrbitGun.hh"
#include "EventGenerator/inc/EjectedNeutronGun.hh"
#include "EventGenerator/inc/EjectedProtonGun.hh"
#include "EventGenerator/inc/EjectedPhotonGun.hh"
#include "EventGenerator/inc/FromG4BLFile.hh"
#include "EventGenerator/inc/NuclearCaptureGun.hh"
#include "EventGenerator/inc/ParticleGun.hh"
#include "EventGenerator/inc/PiCapture.hh"
#include "EventGenerator/inc/PiEplusNuGun.hh"
#include "EventGenerator/inc/PrimaryProtonGun.hh"
#include "EventGenerator/inc/ExtMonFNALGun.hh"
#include "EventGenerator/inc/BremElectronGun.hh"
#include "EventGenerator/inc/MuonicXRayGun.hh"
#include "EventGenerator/inc/CaloCalibGun.hh"
#include "SeedService/inc/SeedService.hh"

// Includes from art and its toolchain.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Other external includes.
#include <boost/shared_ptr.hpp>

// C++ includes.
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

namespace mu2e {

  class EventGenerator : public art::EDProducer {

  public:

    explicit EventGenerator(fhicl::ParameterSet const& pSet);
    // Accept compiler written d'tor.  Modules are never moved or copied.

    virtual void produce (art::Event& e);
    virtual void beginRun(art::Run&   r);

  private:

    // Name of the run-time configuration file.
    string _configfile;

    // Control the behaviour of messages from the SimpleConfig object holding
    // the geometry parameters.
    bool _allowReplacement;
    bool _messageOnReplacement;
    bool _messageOnDefault;
    int  _configStatsVerbosity;

    // Print final config file after all replacements.
    bool _printConfig;

    // The code will only call produces for G4BeamlineInfoCollection
    // when that product will actually be produced.
    bool _produceG4blInfo;

    // A collection of all of the generators that we will run.
    typedef  boost::shared_ptr<GeneratorBase> GeneratorBasePtr;
    std::vector<GeneratorBasePtr> _generators;

    // Check for a valid configuration
    void checkConfig( const SimpleConfig&  config);

    // Check to see if any of the code that will fill the G4BeamlineInfoCollection
    // are present in the configuration.
    bool checkForG4blFile();

  };

  EventGenerator::EventGenerator(fhicl::ParameterSet const& pSet):
    _configfile(           pSet.get<std::string>("inputfile",             "generatorconfig.txt")),
    _allowReplacement(     pSet.get<bool>       ("allowReplacement",      true)),
    _messageOnReplacement( pSet.get<bool>       ("messageOnReplacement",  false)),
    _messageOnDefault(     pSet.get<bool>       ("messageOnDefault",      false)),
    _configStatsVerbosity( pSet.get<int>        ("configStatsVerbosity",  0)),
    _printConfig(          pSet.get<bool>       ("printConfig",           false)),
    _produceG4blInfo(checkForG4blFile()){

    produces<GenParticleCollection>();
    if ( _produceG4blInfo  ){
      produces<G4BeamlineInfoCollection>();
    }

    // A common random engine for the generators to use.
    createEngine( art::ServiceHandle<SeedService>()->getSeed() );
  }

  // Only some of the event generators produce G4BeamlineInfo.
  // Check to see if any are selected.
  bool EventGenerator::checkForG4blFile(){
    SimpleConfig config(_configfile, true, false, false );

    // At the present time the only code that produces G4BeamlineInfo
    // is the FromG4BLFile option.
    return config.getBool( "fromG4BLFile.do", false );
  }

  // Run these at beginRun time in case any of them depend on geoemtry
  // information that may change with conditions information.
  // Otherwise this information could be computed in the c'tor.
  void EventGenerator::beginRun( art::Run &run){

    static int ncalls(0);
    if ( ++ncalls > 1){
      mf::LogInfo("EventGenerator")
        << "EventGenerator does not change state at beginRun.  Hope that's OK.";
      return;
    }

    cout << "Event generator configuration file: "
         << _configfile
         << "\n"
         << endl;

    SimpleConfig config(_configfile, _allowReplacement, _messageOnReplacement, _messageOnDefault );
    checkConfig(config);

    if ( _printConfig ){
      config.print(cout,"EvtGen: ");
    }

    // Change this to modify rather than delete and make an new one??

    // Delete generators from the previous run.
    _generators.clear();

    // Which generators will we run?
    bool doConv                 = config.getBool( "conversionGun.do",    false );
    bool doParticleGun          = config.getBool( "particleGun.do",      false );
    bool doCosmicDYB            = config.getBool( "cosmicDYB.do",        false );
    bool doPiCapture            = config.getBool( "picapture.do",        false );
    bool doEjectedProton        = config.getBool( "ejectedProtonGun.do", false );
    bool doEjectedNeutron       = config.getBool( "ejectedNeutronGun.do",false );
    bool doEjectedPhoton        = config.getBool( "ejectedPhotonGun.do", false );
    bool doDIO                  = config.getBool( "decayinorbitGun.do",  false );
    bool doPiEplusNu            = config.getBool( "piEplusNuGun.do",     false );
    bool doPrimaryProtonGun     = config.getBool( "primaryProtonGun.do", false );
    bool doFromG4BLFile         = config.getBool( "fromG4BLFile.do",     false );
    bool doNuclearCapture       = config.getBool( "nuclearCaptureGun.do",false );
    bool doExtMonFNALGun        = config.getBool( "extMonFNALGun.do",    false );
    bool doStoppedMuonGun       = config.getBool( "stoppedMuonGun.do",   false );
    bool doBremElectronGun      = config.getBool( "bremElectronGun.do",  false );
    bool doMuonicXRayGun        = config.getBool( "muonicXRayGun.do",    false );
    bool doCaloCalibGun         = config.getBool( "caloCalibGun.do",     false );

    // Instantiate generators for this run.
    if ( doParticleGun)          _generators.push_back( GeneratorBasePtr( new ParticleGun(      run, config)) );
    if ( doConv)                 _generators.push_back( GeneratorBasePtr( new ConversionGun(    run, config)) );
    if ( doCosmicDYB)            _generators.push_back( GeneratorBasePtr( new CosmicDYB(        run, config)) );
    if ( doPiCapture)            _generators.push_back( GeneratorBasePtr( new PiCapture(        run, config)) );
    if ( doDIO)                  _generators.push_back( GeneratorBasePtr( new DecayInOrbitGun(  run, config)) );
    if ( doEjectedProton)        _generators.push_back( GeneratorBasePtr( new EjectedProtonGun( run, config)) );
    if ( doEjectedNeutron)       _generators.push_back( GeneratorBasePtr( new EjectedNeutronGun(run, config)) );
    if ( doEjectedPhoton)        _generators.push_back( GeneratorBasePtr( new EjectedPhotonGun (run, config)) );
    if ( doPiEplusNu)            _generators.push_back( GeneratorBasePtr( new PiEplusNuGun(     run, config)) );
    if ( doPrimaryProtonGun)     _generators.push_back( GeneratorBasePtr( new PrimaryProtonGun( run, config)) );
    if ( doFromG4BLFile)         _generators.push_back( GeneratorBasePtr( new FromG4BLFile(     run, config)) );
    if ( doNuclearCapture)       _generators.push_back( GeneratorBasePtr( new NuclearCaptureGun(run, config)) );
    if ( doExtMonFNALGun)        _generators.push_back( GeneratorBasePtr( new ExtMonFNALGun(    run, config)) );
    if ( doStoppedMuonGun)       _generators.push_back( GeneratorBasePtr( new StoppedMuonGun(   run, config)) );
    if ( doBremElectronGun)      _generators.push_back( GeneratorBasePtr( new BremElectronGun(  run, config)) );
    if ( doMuonicXRayGun)        _generators.push_back( GeneratorBasePtr( new MuonicXRayGun(    run, config)) );
    if ( doCaloCalibGun)         _generators.push_back( GeneratorBasePtr( new CaloCalibGun(     run, config)) );

    if ( _generators.size() == 0 ){
      mf::LogWarning("CONTROL")
        << "EventGenerator has no generators enabled. Hope that's OK.";
    }

    config.printAllSummaries( cout, _configStatsVerbosity, "EvtGen: ");

  }

  void
  EventGenerator::produce(art::Event& evt) {

    // Make the collection to hold the output.
    unique_ptr<GenParticleCollection> genParticles(new GenParticleCollection);
    unique_ptr<G4BeamlineInfoCollection> g4beamlineData(new G4BeamlineInfoCollection);

    // Run all of the registered generators.
    for ( std::vector<GeneratorBasePtr>::const_iterator
            i = _generators.begin(),
            e = _generators.end();
          i !=e; ++i ){
            if ((*i)->isFromG4bl()) {
                    ((FromG4BLFile*)i->get())->generate(*genParticles,g4beamlineData.get());
            } else {
                    (*i)->generate(*genParticles);
            }
    }

    // Put the generated particles into the event.
    evt.put(std::move(genParticles));
    if (_produceG4blInfo) {
      evt.put(std::move(g4beamlineData));
    }

  }

  // Look for inconsistencies in the config file.
  void EventGenerator::checkConfig( const SimpleConfig&  config){

    // The known cosmic ray generators.
    vector<string> keys;
    if ( keys.size() == 0 ){
      keys.push_back("cosmictoy.do");
      keys.push_back("cosmicDYB.do");
      keys.push_back("cosmic.do");
    }

    // Require that 0 or 1 of the generators to be present.
    requireUniqueKey( keys, config );

  }


}


using mu2e::EventGenerator;
DEFINE_ART_MODULE(EventGenerator);
