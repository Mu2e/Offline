// A helper class to randomly sample entries read from a collection of
// ROOT trees.  Input ROOT files are read only once, during the
// construction phase, and data are stored in memory.  There are two
// template arguments: EventRecord and NtupleRecord.  They can be the
// same type (the default), then a single row of an input tree will be
// read per EventRecord.   Otherwise EventRecord is expected to be
// a std::vector<NtupleRecord>, and we will use two branches from
// the input trees:  a "main record" branch that corresponds to
// NtupleRecord, and a "particle in event (pie)" branch that
// stores a single integer, and allows to keep correlations
// between particles (NtupleRecord) in an event (EventRecord).
//
// The fire() method returns one of the stored EventRecord entries.
// All the stored entries are sampled with equal probabilities.
//
// The optional averageNumRecordsToUse parameter controls memory use:
// if the number of input records exceeds the parameter, not all
// records are stored.  Entries to be used are randomly selected at
// the initialization time with equal *per-event* probabilities to
// make the expected number of stored ntuple records *approximately*
// equal to the specified parameter.  The approximations should be
// good if the average per-event ntuple record multiplicity is not
// large.  The biases affect how much memory the code uses, but do not
// affect physics results.  Setting averageNumRecordsToUse<=0 disables
// the feature (the default). If the number of inputs is less than
// averageNumRecordsToUse, all input records are used.
//
// See StoppedParticleReactionGun_module.cc and InFlightParticleSampler_module.cc
// for examples of use.
//
// Andrei Gaponenko, 2014, 2015

#ifndef RootTreeSampler_hh
#define RootTreeSampler_hh

#include <vector>

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/TupleAs.h"
#include "fhiclcpp/ParameterSet.h"

#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandFlat.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"

#include "cetlib_except/exception.h"

#include "TTree.h"
#include "TFile.h"

#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"

namespace mu2e {

  template<class EventRecord, class NtupleRecord=EventRecord>
  class RootTreeSampler {
  public:

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Sequence<std::string> inputFiles {
        Name("inputFiles"),
          Comment("List of input ntuple files")
      };

      fhicl::Atom<std::string> treeName {
        Name("treeName"),
          Comment("Name of the ROOT tree object containing particle records")
          };

      fhicl::Atom<std::string> branchName {
        Name("branchName"),
          Comment("Name of the ROOT tree branch containing particle records")
          };

      fhicl::Atom<long> averageNumRecordsToUse {
        Name("averageNumRecordsToUse"),
          Comment("Zero means use all record.  Use a non-zero value\n"
                  "to load only a random subset of the records to limit memory use."),
          0
          };

      fhicl::Atom<long> verbosityLevel {
        Name("verbosityLevel"),
          Comment("A positive value generates more printouts than the default."),
          0
          };

      fhicl::Atom<std::string> pieBranchName {
        Name("pieBranchName"),
          Comment("Only for generators resampling correlated particles from\n"
                  "multiple ntuple records:\n"
                  "Name of the ROOT tree branch that maps input records to events."),
          [this](){ return !std::is_same<EventRecord,NtupleRecord>::value; }
      };
    };

    RootTreeSampler(art::RandomNumberGenerator::base_engine_t& engine,
                    const Config& conf);

    // legacy interface with pset
    RootTreeSampler(art::RandomNumberGenerator::base_engine_t& engine,
                    const fhicl::ParameterSet& pset);

    const EventRecord& fire() { return records_.at(randFlat_.fireInt(records_.size())); }

    typename std::vector<EventRecord>::size_type
    numRecords() const { return records_.size(); }

  private:
    CLHEP::RandFlat randFlat_;
    std::vector<EventRecord> records_;

    typedef std::vector<std::string> Strings;

    long countInputRecords(const art::ServiceHandle<art::TFileService>& tfs,
                           const Strings& files,
                           const std::string& treeName);

    //----------------------------------------------------------------
    struct SingleRecordGetter {
      RootTreeSampler *parent_;
      TTree *nt_;
      TBranch *mrb_;
      NtupleRecord ntr_;

      long averageNumRecordsToUse_;
      double recordUseFraction_;

      bool hasMore_ = false;
      Long64_t numTreeEntries_;
      Long64_t currentEntry_ = 0;
      Long64_t numUsedNtupleEntries_ = 0;

      SingleRecordGetter(RootTreeSampler *ts,
                         TTree *tr,
                         TBranch *mainRecordBranch,
                         const Config& conf,
                         double recordUseFraction)
        : parent_(ts)
        , nt_(tr)
        , mrb_(mainRecordBranch)
        , recordUseFraction_(recordUseFraction)
        , numTreeEntries_(tr->GetEntries())
      {
        mrb_->SetAddress(&ntr_);
      }

      SingleRecordGetter(RootTreeSampler *ts,
                         TTree *tr,
                         TBranch *mainRecordBranch,
                         const fhicl::ParameterSet& pset,
                         double recordUseFraction)
        : parent_(ts)
        , nt_(tr)
        , mrb_(mainRecordBranch)
        , recordUseFraction_(recordUseFraction)
        , numTreeEntries_(tr->GetEntries())
      {
        mrb_->SetAddress(&ntr_);
      }

      bool hasMoreRecords() {
        while(currentEntry_ < numTreeEntries_) {
          if(parent_->randFlat_.fire() < recordUseFraction_) {
            mrb_->GetEntry(currentEntry_++);
            ++numUsedNtupleEntries_;
            return hasMore_ = true;
          }
          ++currentEntry_;
        }
        return hasMore_ = false;
      }

      EventRecord getRecord() const {
        if(!hasMore_) {
          throw cet::exception("BUG")
            <<"RootTreeSampler::SingleRecordGetter(): getRecord() called with no current record.\n";

        }
        return ntr_;
      }

      Long64_t numUsedNtupleEntries() const { return numUsedNtupleEntries_; }

    }; // SingleRecordGetter

    //----------------------------------------------------------------
    struct MultiRecordGetter {
      RootTreeSampler *parent_;
      TTree *nt_;
      TBranch *mrb_;
      NtupleRecord ntr_;
      EventRecord evt_;

      TBranch *pieb_;
      unsigned particleInEvent_;

      long averageNumRecordsToUse_;
      double recordUseFraction_;

      Long64_t numTreeEntries_;
      Long64_t currentEntry_;

      Long64_t numUsedNtupleEntries_;

      MultiRecordGetter(RootTreeSampler *ts,
                        TTree *tr,
                        TBranch *mainRecordBranch,
                        const Config& conf,
                        double recordUseFraction)
        : parent_(ts)
        , nt_(tr)
        , mrb_(mainRecordBranch)

        , pieb_(nullptr)
        , particleInEvent_(-1)

        , recordUseFraction_(recordUseFraction)
        , numTreeEntries_(tr->GetEntries())
        , currentEntry_(0)
        , numUsedNtupleEntries_(0)
      {
        const std::string pbn = conf.pieBranchName();
        pieb_ = tr->GetBranch(pbn.c_str());

        mrb_->SetAddress(&ntr_);

        if(!pieb_) {
          throw cet::exception("BADINPUT")
            <<"RootTreeSampler: Could not get branch \""<<conf.pieBranchName()
            <<"\" in tree \""<<nt_->GetName()
            <<"\"\n";
        }
        pieb_->SetAddress(&particleInEvent_);
      }

      MultiRecordGetter(RootTreeSampler *ts,
                        TTree *tr,
                        TBranch *mainRecordBranch,
                        const fhicl::ParameterSet& pset,
                        double recordUseFraction)
        : parent_(ts)
        , nt_(tr)
        , mrb_(mainRecordBranch)

        , pieb_(nullptr)
        , particleInEvent_(-1)

        , recordUseFraction_(recordUseFraction)
        , numTreeEntries_(tr->GetEntries())
        , currentEntry_(0)
        , numUsedNtupleEntries_(0)
      {
        const std::string pbn = pset.get<std::string>("pieBranchName");
        pieb_ = tr->GetBranch(pbn.c_str());

        mrb_->SetAddress(&ntr_);

        if(!pieb_) {
          throw cet::exception("BADINPUT")
            <<"RootTreeSampler: Could not get branch \""<<pset.get<std::string>("pieBranchName")
            <<"\" in tree \""<<nt_->GetName()
            <<"\"\n";
        }
        pieb_->SetAddress(&particleInEvent_);
      }

      bool hasMoreRecords() {
        evt_.clear();

        while(currentEntry_ < numTreeEntries_) {
          pieb_->GetEntry(currentEntry_);
          if(particleInEvent_ != 0) { // something went wrong - we should be aligned on the beginning of an event here
            throw cet::exception("BADINPUT")<<"RootTreeSampler: Error: unexpected particleInEvent!=0";
          }

          const bool use = (parent_->randFlat_.fire() < recordUseFraction_);
          if(use) { // read the current record
            mrb_->GetEntry(currentEntry_);
            evt_.emplace_back(ntr_);
            ++numUsedNtupleEntries_;
          }

          // go through the rest of records from this event
          for(++currentEntry_; currentEntry_ < numTreeEntries_; ++currentEntry_) {

            pieb_->GetEntry(currentEntry_);

            if(particleInEvent_ == 0) {
              // that's the first record from the next even.
              // the current event is complete
              if(use) {
                return true;
              }
              else {
                break; // We did not use that event.  Go to the next one in the outer loop.
              }
            }

            if(use) {
              mrb_->GetEntry(currentEntry_);
              evt_.emplace_back(ntr_);
              ++numUsedNtupleEntries_;
            }
          }
        }

        return false;
      }

      EventRecord getRecord() const {
        if(evt_.empty()) {
          throw cet::exception("BUG")
            <<"RootTreeSampler::SingleRecordGetter(): "
            <<"getRecord() called with no current record.\n";
        }
        return evt_;
      }

      Long64_t numUsedNtupleEntries() const { return numUsedNtupleEntries_; }

    }; // MultiRecordGetter
    //----------------------------------------------------------------

  }; // RootTreeSampler
}

//================================================================
namespace mu2e {

  template<class EventRecord, class NtupleRecord>
  RootTreeSampler<EventRecord, NtupleRecord>::
  RootTreeSampler(art::RandomNumberGenerator::base_engine_t& engine,
                  const Config& conf)
    : randFlat_(engine)
  {
    const auto inputFiles(conf.inputFiles());
    const auto treeName(conf.treeName());
    const long averageNumRecordsToUse(conf.averageNumRecordsToUse());
    int verbosityLevel = conf.verbosityLevel();

    if(inputFiles.empty()) {
      throw cet::exception("BADCONFIG")<<"Error: no inputFiles";
    }

    art::ServiceHandle<art::TFileService> tfs;

    double recordUseFraction = 1.;
    if(averageNumRecordsToUse > 0) {
      const long totalRecords = countInputRecords(tfs, inputFiles, treeName);
      if(averageNumRecordsToUse < totalRecords) {
        recordUseFraction = averageNumRecordsToUse/double(totalRecords);
        if(verbosityLevel > 0) {
          std::cout<<"RootTreeSampler: recordUseFraction = "<<recordUseFraction
                   <<" for the requested average = "<<averageNumRecordsToUse
                   <<" and total number of input records = "<<totalRecords
                   <<std::endl;
        }
      }
      else {
        if(verbosityLevel > 0) {
          std::cout<<"RootTreeSampler: forcing recordUseFraction=1 "
                   <<": the requested number = "<<averageNumRecordsToUse
                   <<" exceeds the number of available inputs = "<<totalRecords
                   <<std::endl;
        }
      }
    }

    //----------------------------------------------------------------
    // Load the records
    for(const auto& fn : inputFiles) {
      const std::string resolvedFileName = ConfigFileLookupPolicy()(fn);
      TFile *infile = tfs->make<TFile>(resolvedFileName.c_str(), "READ");

      TTree *nt = dynamic_cast<TTree*>(infile->Get(treeName.c_str()));
      if(!nt) {
        throw cet::exception("BADINPUT")<<"RootTreeSampler: Could not get tree \""<<treeName
                                        <<"\" from file \""<<infile->GetName()
                                        <<"\"\n";
      }

      //----------------------------------------------------------------
      // A rudimentary check that the input ntuple is consistent with the data structure

      const auto branchName(conf.branchName());
      TBranch *bb = nt->GetBranch(branchName.c_str());
      if(!bb) {
        throw cet::exception("BADINPUT")<<"RootTreeSampler: Could not get branch \""<<branchName
                                        <<"\" in tree \""<<treeName
                                        <<"\" from file \""<<infile->GetName()
                                        <<"\"\n";
      }

      if(unsigned(bb->GetNleaves()) != NtupleRecord::numBranchLeaves()) {
        throw cet::exception("BADINPUT")<<"RootTreeSampler: wrong number of leaves: expect "
                                        <<NtupleRecord::numBranchLeaves()<<", but branch \""<<branchName
                                        <<"\", tree \""<<treeName
                                        <<"\" in file \""<<infile->GetName()
                                        <<"\" has "<<bb->GetNleaves()
                                        <<"\n";
      }
      //----------------------------------------------------------------

      const Long64_t nTreeEntries = nt->GetEntries();
      std::cout<<"RootTreeSampler: reading "<<nTreeEntries
               <<" entries.  Tree "<<treeName
               <<", file "<<infile->GetName()
               <<std::endl;

      // If the average per-event ntuple record multiplicity is large,
      // this will over-allocate the memory.
      records_.reserve(records_.size()
                       // Add "mean + 3sigma": do not re-allocate in most cases.
                       + recordUseFraction*nTreeEntries
                       + 3*recordUseFraction*sqrt(double(nTreeEntries)));

      typename std::conditional
        <std::is_same<EventRecord,NtupleRecord>::value,
        SingleRecordGetter, MultiRecordGetter>::type
        egt(this, nt, bb, conf, recordUseFraction);

      while(egt.hasMoreRecords()) {
        records_.push_back(egt.getRecord());
      }

      if(verbosityLevel > 0) {
        std::cout<<"RootTreeSampler: stored "<<records_.size()
                 <<" event entries.  Used "<<egt.numUsedNtupleEntries()
                 <<" ntuple entries."
                 <<std::endl;
      }

    } // for(inputFiles)

  } // Constructor (conf)

  //================================================================
  template<class EventRecord, class NtupleRecord>
  RootTreeSampler<EventRecord, NtupleRecord>::
  RootTreeSampler(art::RandomNumberGenerator::base_engine_t& engine,
                  const fhicl::ParameterSet& pset)
    : randFlat_(engine)
  {
    const auto inputFiles(pset.get<std::vector<std::string> >("inputFiles"));
    const auto treeName(pset.get<std::string>("treeName"));
    const long averageNumRecordsToUse(pset.get<long>("averageNumRecordsToUse", 0));
    int verbosityLevel = pset.get<int>("verbosityLevel", 0);

    if(inputFiles.empty()) {
      throw cet::exception("BADCONFIG")<<"Error: no inputFiles";
    }

    art::ServiceHandle<art::TFileService> tfs;

    double recordUseFraction = 1.;
    if(averageNumRecordsToUse > 0) {
      const long totalRecords = countInputRecords(tfs, inputFiles, treeName);
      if(averageNumRecordsToUse < totalRecords) {
        recordUseFraction = averageNumRecordsToUse/double(totalRecords);
        if(verbosityLevel > 0) {
          std::cout<<"RootTreeSampler: recordUseFraction = "<<recordUseFraction
                   <<" for the requested average = "<<averageNumRecordsToUse
                   <<" and total number of input records = "<<totalRecords
                   <<std::endl;
        }
      }
      else {
        if(verbosityLevel > 0) {
          std::cout<<"RootTreeSampler: forcing recordUseFraction=1 "
                   <<": the requested number = "<<averageNumRecordsToUse
                   <<" exceeds the number of available inputs = "<<totalRecords
                   <<std::endl;
        }
      }
    }

    //----------------------------------------------------------------
    // Load the records
    for(const auto& fn : inputFiles) {
      const std::string resolvedFileName = ConfigFileLookupPolicy()(fn);
      TFile *infile = tfs->make<TFile>(resolvedFileName.c_str(), "READ");

      TTree *nt = dynamic_cast<TTree*>(infile->Get(treeName.c_str()));
      if(!nt) {
        throw cet::exception("BADINPUT")<<"RootTreeSampler: Could not get tree \""<<treeName
                                        <<"\" from file \""<<infile->GetName()
                                        <<"\"\n";
      }

      //----------------------------------------------------------------
      // A rudimentary check that the input ntuple is consistent with the data structure

      const auto branchName(pset.get<std::string>("branchName"));
      TBranch *bb = nt->GetBranch(branchName.c_str());
      if(!bb) {
        throw cet::exception("BADINPUT")<<"RootTreeSampler: Could not get branch \""<<branchName
                                        <<"\" in tree \""<<treeName
                                        <<"\" from file \""<<infile->GetName()
                                        <<"\"\n";
      }

      if(unsigned(bb->GetNleaves()) != NtupleRecord::numBranchLeaves()) {
        throw cet::exception("BADINPUT")<<"RootTreeSampler: wrong number of leaves: expect "
                                        <<NtupleRecord::numBranchLeaves()<<", but branch \""<<branchName
                                        <<"\", tree \""<<treeName
                                        <<"\" in file \""<<infile->GetName()
                                        <<"\" has "<<bb->GetNleaves()
                                        <<"\n";
      }
      //----------------------------------------------------------------

      const Long64_t nTreeEntries = nt->GetEntries();
      std::cout<<"RootTreeSampler: reading "<<nTreeEntries
               <<" entries.  Tree "<<treeName
               <<", file "<<infile->GetName()
               <<std::endl;

      // If the average per-event ntuple record multiplicity is large,
      // this will over-allocate the memory.
      records_.reserve(records_.size()
                       // Add "mean + 3sigma": do not re-allocate in most cases.
                       + recordUseFraction*nTreeEntries
                       + 3*recordUseFraction*sqrt(double(nTreeEntries)));

      typename std::conditional
        <std::is_same<EventRecord,NtupleRecord>::value,
        SingleRecordGetter, MultiRecordGetter>::type
        egt(this, nt, bb, pset, recordUseFraction);

      while(egt.hasMoreRecords()) {
        records_.push_back(egt.getRecord());
      }

      if(verbosityLevel > 0) {
        std::cout<<"RootTreeSampler: stored "<<records_.size()
                 <<" event entries.  Used "<<egt.numUsedNtupleEntries()
                 <<" ntuple entries."
                 <<std::endl;
      }

    } // for(inputFiles)

  } // Constructor (pset)

  //================================================================
  template<class EventRecord, class NtupleRecord>
  long RootTreeSampler<EventRecord, NtupleRecord>::countInputRecords(const art::ServiceHandle<art::TFileService>& tfs,
                                                                     const Strings& inputFiles,
                                                                     const std::string& treeName)
  {
    long res=0;
    for(const auto& fn : inputFiles) {
      const std::string resolvedFileName = ConfigFileLookupPolicy()(fn);
      TFile *infile = tfs->make<TFile>(resolvedFileName.c_str(), "READ");
      TTree *nt = dynamic_cast<TTree*>(infile->Get(treeName.c_str()));
      if(!nt) {
        throw cet::exception("BADINPUT")<<"Could not get tree \""<<treeName
                                        <<"\" from file \""<<infile->GetName()
                                        <<"\"\n";
      }

      const Long64_t nTreeEntries = nt->GetEntries();
      res += nTreeEntries;
    }
    return res;
  }

  //================================================================
}

#endif /* RootTreeSampler_hh */
