//
// This source reads a text-formated sequence of #s of protons at the production target
// and translate them into a set of events containing the corresponding PBI number.
// Original author: David Bown (LBNL), 2020
//
// The SubRun product of type ProtonBunchIntensity has been deprecated; it is retained
// for compatibilty with the existing MDC2020DEV files and will be removed at a convenient
// time.  It's successor is the Subrun product ProtonBunchIntensitySummary.
//
// If the parameter fileNames has more than one file then
//   - a new subrun will be started for each file
//   - the subrun number will be incremented by 1 for each file
//   - the code is configurable so that the subrun data products can be computed
//     separately for each subrun or integrated over all subruns.
//
// There is an option, integratedSummary that will
//   - read all input files at the start of the job and compute the
//     subrun summary information using all input files.
//   - this will be added to the subrun scope data products and will
//     be the same in all subruns.
//
// Fixme:
//   Are we treating the RunPrincipal correctly?
//     - Suppose we create a run product at the start of the first subrun.
//     - When we start a new subrun in the same run, I think that the run
//       product disappears because we made a new RunPrincipal.
//   So long as we do not write any run products this is moot.
//   So long as we do not write multiple subruns to same output file it is also moot.
//

#include <iostream>
#include <fstream>
#include <boost/utility.hpp>
#include <cassert>
#include <set>
#include <string>
#include <cstdio>
#include <vector>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>

#include "fhiclcpp/types/Name.h"
#include "art/Framework/IO/Sources/Source.h"
#include "art/Framework/Core/InputSourceMacros.h"
#include "art/Framework/IO/Sources/SourceHelper.h"
#include "art/Framework/Principal/RunPrincipal.h"
#include "art/Framework/Principal/SubRunPrincipal.h"
#include "art/Framework/Principal/EventPrincipal.h"
#include "art/Framework/IO/Sources/put_product_in_principal.h"
#include "canvas/Persistency/Provenance/Timestamp.h"
#include "canvas/Persistency/Provenance/RunID.h"
#include "canvas/Persistency/Provenance/SubRunID.h"
#include "canvas/Persistency/Provenance/EventID.h"
#include "canvas/Persistency/Provenance/BranchType.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "MCDataProducts/inc/ProtonBunchIntensity.hh"
#include "MCDataProducts/inc/ProtonBunchIntensitySummary.hh"
using namespace std;

namespace mu2e {
  using namespace boost::accumulators;

  struct Config
  {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    fhicl::Sequence<std::string> inputFiles{Name("fileNames"),Comment("List of text-formated PBI sequence files")};
    fhicl::Atom<std::string> moduleLabel{Name("reconstitutedModuleLabel"),Comment("Module label for the reconstituted data products"),"PBISequence"};
    fhicl::Atom<std::string> instance{Name("deprecatedInstanceName"),Comment("Instance name for deprecated subrun data product ProtonBunchIntensity"),""};
    fhicl::Atom<unsigned int> runNumber{Name("runNumber"), Comment("First run number")};
    fhicl::Atom<bool>  integratedSummary{Name("integratedSummary"), Comment("If true, then PBI summary is summed over all input files; else per input file.")};
    fhicl::Atom<unsigned int> verbosity{Name("verbosity"), Comment("Verbosity level; larger means more printout"), 0};

    // These are used by art and are required.
    fhicl::Atom<std::string> module_label{Name("module_label"), Comment("Art module label"), ""};
    fhicl::Atom<std::string> module_type{Name("module_type"), Comment("Art module type"), ""};
  };
  typedef fhicl::WrappedTable<Config> Parameters;

 //================================================================
  class PBISequenceDetail : private boost::noncopyable {
    std::string myModuleLabel_;
    std::string deprecatedInstanceName_;
    art::SourceHelper const& pm_;
    unsigned runNumber_;
    unsigned subRunNumber_;
    unsigned currentEventNumber_;
    unsigned verbosity_;
    art::SubRunID lastSubRunID_;
    accumulator_set<double, stats< tag::mean, tag::variance>> nprotAcc_; // accumulator for proton statistics for one input file (subrun).
    accumulator_set<double, stats< tag::mean, tag::variance>> runNprotAcc_; // accumulator for proton statistics over the full input file set (run).

    std::string currentFileName_;
    std::ifstream *currentFile_ = nullptr;

    bool  integratedSummary_;

    // A helper function used to manage the principals.
    // This is boilerplate that does not change if you change the data products.
    void managePrincipals ( int runNumber,
	int subRunNumber,
	int eventNumber,
	art::RunPrincipal*&    outR,
	art::SubRunPrincipal*& outSR,
	art::EventPrincipal*&  outE);

    void computeRunDataProducts( std::vector<std::string> const& inputFiles );

    public:
    PBISequenceDetail(const Parameters &conf,
	art::ProductRegistryHelper &,
	const art::SourceHelper &);

      void readFile(std::string const& filename, art::FileBlock*& fb);

      bool readNext(art::RunPrincipal* const& inR,
                    art::SubRunPrincipal* const& inSR,
                    art::RunPrincipal*& outR,
                    art::SubRunPrincipal*& outSR,
                    art::EventPrincipal*& outE);

      void closeCurrentFile();

    };

    //----------------------------------------------------------------
    PBISequenceDetail::PBISequenceDetail(const Parameters& conf,
                                             art::ProductRegistryHelper& rh,
                                             const art::SourceHelper& pm)
      : myModuleLabel_(conf().moduleLabel())
      , deprecatedInstanceName_(conf().instance())
      , pm_(pm)
      , runNumber_(conf().runNumber())
      , subRunNumber_(-1U)
      , currentEventNumber_(0)
      , verbosity_(conf().verbosity())
      , integratedSummary_(conf().integratedSummary())
    {
      if(!art::RunID(runNumber_).isValid()) {
        throw cet::exception("BADCONFIG", " PBISequence: ")
          << " fhicl::ParameterSet specifies an invalid runNumber = "<<runNumber_<<"\n";
      }
//    For a source module, "produces" is spelled "reconstitutes"
      rh.reconstitutes<mu2e::ProtonBunchIntensity,art::InEvent>(myModuleLabel_);
      rh.reconstitutes<mu2e::ProtonBunchIntensity,art::InSubRun>(myModuleLabel_,deprecatedInstanceName_);
      rh.reconstitutes<mu2e::ProtonBunchIntensitySummary,art::InSubRun>(myModuleLabel_);

      computeRunDataProducts( conf().inputFiles() );
    }

    //----------------------------------------------------------------
    void PBISequenceDetail::readFile(const std::string& filename, art::FileBlock*& fb) {
      using namespace boost::accumulators;
      currentFileName_ = filename;
      currentFile_ = new ifstream(currentFileName_,std::ifstream::in);
      // reset counters
      subRunNumber_++;
      currentEventNumber_ = 0;
      // compute statistics on protons in this file
      nprotAcc_ = {};
      double protons;
      unsigned nprotons;
      while(true) {
	*currentFile_ >> protons;
	if ( currentFile_->good()) {
	  nprotons = (unsigned)rint(protons);
	  nprotAcc_(double(nprotons));
	} else
	break;
      }
      if ( verbosity_ > 0 ){
        std::cout << "Read " << extract_result<tag::count>(nprotAcc_) << " events with " << extract_result<tag::mean>(nprotAcc_) << " <protons> " << extract_result<tag::variance>(nprotAcc_) << " variance from file " << currentFileName_ << std::endl;
      }
      // rewind the file
      currentFile_->clear(ios::eofbit);
      currentFile_->seekg (0, ios::beg);
      fb = new art::FileBlock(art::FileFormatVersion(1, "PBISequenceTextInput"), currentFileName_);
    }

    //----------------------------------------------------------------
    void PBISequenceDetail::closeCurrentFile() {
      currentFileName_ = "";
      currentFile_->close();
      delete currentFile_;
      currentFile_ = nullptr;
    }

    //----------------------------------------------------------------
    bool PBISequenceDetail::readNext(art::RunPrincipal* const& inR,
                                       art::SubRunPrincipal* const& inSR,
                                       art::RunPrincipal*& outR,
                                       art::SubRunPrincipal*& outSR,
                                       art::EventPrincipal*& outE)
    {
      double protons;
      unsigned nprotons;
      (*currentFile_) >>  protons;
      if (!currentFile_->good()) return false;
      nprotons = (unsigned)rint(protons);
      managePrincipals(runNumber_, subRunNumber_, ++currentEventNumber_, outR, outSR, outE);
      std::unique_ptr<ProtonBunchIntensity> pbi(new ProtonBunchIntensity(nprotons));
      art::put_product_in_principal(std::move(pbi), *outE, myModuleLabel_);
      return true;

    } // readNext()


  // Each time that we encounter a new run, a new subRun or a new event, we need to make a new principal
  // of the appropriate type.  This code does not need to change as the number and type of data products changes.
  void PBISequenceDetail::managePrincipals ( int runNumber,
                                          int subRunNumber,
                                          int eventNumber,
                                          art::RunPrincipal*&    outR,
                                          art::SubRunPrincipal*& outSR,
                                          art::EventPrincipal*&  outE){

    art::Timestamp ts;

    art::SubRunID newID(runNumber, subRunNumber);

    if(newID != lastSubRunID_) {
      outR = pm_.makeRunPrincipal(runNumber, ts);
      // art takes ownership of the object pointed to by outSR and will delete it at the appropriate time.
      outSR = pm_.makeSubRunPrincipal(runNumber,
                                      subRunNumber,
                                      ts);
      // create the PBI subrun summary object
      std::unique_ptr<ProtonBunchIntensitySummary> pbis = ( integratedSummary_ ) ?
        std::make_unique<ProtonBunchIntensitySummary>( extract_result<tag::count>(runNprotAcc_),extract_result<tag::mean>(runNprotAcc_), extract_result<tag::variance>(runNprotAcc_)) :
        std::make_unique<ProtonBunchIntensitySummary>( extract_result<tag::count>(nprotAcc_),extract_result<tag::mean>(nprotAcc_), extract_result<tag::variance>(nprotAcc_));

      if ( verbosity_ > 0 ) {
        std::cout << "SubRun " << subRunNumber << " PBI has SDF = " << pbis->spillDutyFactor() << std::endl;
      }
      art::put_product_in_principal(std::move(pbis), *outSR, myModuleLabel_);
      // for backwards-compatibility
      std::unique_ptr<ProtonBunchIntensity> pbi = ( integratedSummary_ ) ?
        std::make_unique<ProtonBunchIntensity>((unsigned)rint(extract_result<tag::mean>(runNprotAcc_))):
        std::make_unique<ProtonBunchIntensity>((unsigned)rint(extract_result<tag::mean>(nprotAcc_)));
      art::put_product_in_principal(std::move(pbi), *outSR, myModuleLabel_, deprecatedInstanceName_);

    }
    lastSubRunID_ = newID;

    // art takes ownership of the object pointed to by outE and will delete it at the appropriate time.
    outE = pm_.makeEventPrincipal(runNumber, subRunNumber, eventNumber, ts, false);

  } // managePrincipals()


  //----------------------------------------------------------------
  void PBISequenceDetail::computeRunDataProducts( std::vector<std::string> const& inputFiles ){
    for ( auto const& file : inputFiles){
      std::ifstream in(file,std::ifstream::in);
      double protons;
      int nprotons;
      while (in){
	in >> protons;
	if ( in.good()) {
	  nprotons = (unsigned)rint(protons);
	  runNprotAcc_(double(nprotons));
	} else {
          break;
        }
      }
    }
  } // end computeRunDataProducts

} // namespace mu2e

DEFINE_ART_INPUT_SOURCE(art::Source<mu2e::PBISequenceDetail>)
