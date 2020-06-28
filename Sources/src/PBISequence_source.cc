// 
// This source reads a text-formated sequence of #s of protons at the production target
// and translate them into a set of events containing the corresponding PBI number.
// Original author: David Bown (LBNL), 2020

#include <iostream>
#include <fstream>
#include <boost/utility.hpp>
#include <cassert>
#include <set>
#include <string>
#include <cstdio>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>

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
    fhicl::Atom<unsigned int> runNumber{Name("runNumber"), Comment("First run number"), 0};
    fhicl::Atom<std::string> module_label{Name("module_label"), Comment("Art module label"), ""};
    fhicl::Atom<std::string> module_type{Name("module_type"), Comment("Art module type"), ""};
  };
  typedef fhicl::WrappedTable<Config> Parameters;

  //================================================================
  class PBISequenceDetail : private boost::noncopyable {
    std::string myModuleLabel_;
    art::SourceHelper const& pm_;
    unsigned runNumber_;
    unsigned subRunNumber_;
    unsigned currentEventNumber_;
    art::SubRunID lastSubRunID_;
    accumulator_set<float, stats< tag::mean, tag::variance>> nprotAcc_; // accumulator for proton statistics

    std::string currentFileName_;
    std::ifstream *currentFile_ = nullptr;

    // A helper function used to manage the principals.
    // This is boilerplate that does not change if you change the data products.
    void managePrincipals ( int runNumber,
	int subRunNumber,
	int eventNumber,
	art::RunPrincipal*&    outR,
	art::SubRunPrincipal*& outSR,
	art::EventPrincipal*&  outE);


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
      : myModuleLabel_("PBISequence")
      , pm_(pm)
      , runNumber_(conf().runNumber())
      , subRunNumber_(-1U)
      , currentEventNumber_(0)
    {
      if(!art::RunID(runNumber_).isValid()) {
        throw cet::exception("BADCONFIG", " PBISequence: ")
          << " fhicl::ParameterSet specifies an invalid runNumber = "<<runNumber_<<"\n";
      }
// DNB
// I don't understand why 'reconstitutes' is needed instead of 'produces', but it's
// an emperical fact that this module fails in a cryptic way using 'produces'
//
      rh.reconstitutes<mu2e::ProtonBunchIntensity,art::InEvent>(myModuleLabel_);
      rh.reconstitutes<mu2e::ProtonBunchIntensity,art::InSubRun>(myModuleLabel_);
      rh.reconstitutes<mu2e::ProtonBunchIntensitySummary,art::InSubRun>(myModuleLabel_);
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
      float protons;
      unsigned nprotons;
      while(true) {
	*currentFile_ >> protons;
	if ( currentFile_->good()) {
	  nprotons = (unsigned)rint(protons);
	  nprotAcc_(float(nprotons));
	} else
	break;
      }
      std::cout << "Read " << extract_result<tag::count>(nprotAcc_) << " events with " << extract_result<tag::mean>(nprotAcc_) << " <protons> " << extract_result<tag::variance>(nprotAcc_) << " variance from file " << currentFileName_ << std::endl;
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
      float protons;
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
      // art takes ownership of the object pointed to by outSR and will delete it at the appropriate time.
      outR = pm_.makeRunPrincipal(runNumber, ts);
      outSR = pm_.makeSubRunPrincipal(runNumber,
                                      subRunNumber,
                                      ts);
      // create the PBI subrun summary object
      std::unique_ptr<ProtonBunchIntensitySummary> pbis(new ProtonBunchIntensitySummary(extract_result<tag::count>(nprotAcc_),extract_result<tag::mean>(nprotAcc_), extract_result<tag::variance>(nprotAcc_)));
      std::cout << "SubRun " << subRunNumber << " PBI has SDF = " << pbis->spillDutyFactor() << std::endl;
      art::put_product_in_principal(std::move(pbis), *outSR, myModuleLabel_);
      // for backwards-compatibility
      std::unique_ptr<ProtonBunchIntensity> pbi(new ProtonBunchIntensity((unsigned)rint(extract_result<tag::mean>(nprotAcc_))));
      art::put_product_in_principal(std::move(pbi), *outSR, myModuleLabel_);

    }
    lastSubRunID_ = newID;

    // art takes ownership of the object pointed to by outE and will delete it at the appropriate time.
    outE = pm_.makeEventPrincipal(runNumber, subRunNumber, eventNumber, ts, false);

  } // managePrincipals()
    //----------------------------------------------------------------

} // namespace mu2e

typedef art::Source<mu2e::PBISequenceDetail> PBISequence;
DEFINE_ART_INPUT_SOURCE(PBISequence);
