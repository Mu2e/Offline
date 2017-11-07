//
// Maintain up to date Alignment information and serve it to
// other services and to the modules.
//
// Original author David Norvil Brown (UofL), based on
// GeometryService by Rob Kutschke
//

// C++ include files
#include <iostream>
#include <utility>

// Framework include files
#include "canvas/Persistency/Provenance/ModuleDescription.h"
#include "canvas/Persistency/Provenance/EventID.h"
#include "canvas/Persistency/Provenance/Timestamp.h"
#include "canvas/Persistency/Provenance/SubRunID.h"
#include "canvas/Persistency/Provenance/RunID.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e include files
#include "Alignment/inc/AlignmentService.hh"
#include "Alignment/inc/AlignmentMap.hh"

//#include "GeometryService/inc/DetectorSystem.hh"
//#include "GeometryService/src/DetectorSystemMaker.hh"

using namespace std;

namespace mu2e {


  AlignmentService::AlignmentService(fhicl::ParameterSet const& pset,
                                   art::ActivityRegistry&iRegistry) :
    _inputfile(            pset.get<std::string> ("inputFile",            "alignment000.txt")),
    _allowReplacement(     pset.get<bool>        ("allowReplacement",     true)),
    _messageOnReplacement( pset.get<bool>        ("messageOnReplacement", false)),
    _messageOnDefault(     pset.get<bool>        ("messageOnDefault",     false)),
    _configStatsVerbosity( pset.get<int>         ("configStatsVerbosity", 0)),
    _printConfig(          pset.get<bool>        ("printConfig",          false)),
    _config(nullptr),
    _alignmentMap(new AlignmentMap()),
    _run_count()
  {
    iRegistry.sPreBeginRun.watch(this, &AlignmentService::preBeginRun);
    iRegistry.sPostEndJob.watch (this, &AlignmentService::postEndJob );
  }


  void
  AlignmentService::preBeginRun(art::Run const &) {

    if(++_run_count > 1) {
      mf::LogWarning("GEOM") << "This test version does not change alignment on run boundaries.";
      return;
    }

    cout  << "Alignment input file is: " << _inputfile << "\n";

    _config = unique_ptr<SimpleConfig>(new SimpleConfig(_inputfile,
                                                      _allowReplacement,
                                                      _messageOnReplacement,
                                                      _messageOnDefault ));

    // Print final state of file after all substitutions.
    if ( _printConfig      ){ _config->print(cout, "Alignment: ");       }


    if (_config->getBool("hasAlignment",false)) {
      //      _alignmentMap = unique_ptr<AlignmentMap>( new AlignmentMap() );
      _alignmentMap->make( *_config );
    }


  } // preBeginRun()


  // Called after all modules have completed their end of job.
  void   AlignmentService::postEndJob(){
    _config->printAllSummaries( cout, _configStatsVerbosity, "Align: " );
  }



} // end namespace mu2e

DEFINE_ART_SERVICE(mu2e::AlignmentService);
