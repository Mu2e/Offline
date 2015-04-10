// A helper class to randomly sample entries read from a collection of
// ROOT trees.  It is templated on a structure that must correspond to
// a tree branch to be read.  Input ROOT files are read only once,
// during the construction phase, and data are stored in memory.  The
// optional averageNumRecordsToUse parameter controls memory use: if
// the number of input records exceeds the parameter, not all records
// are stored.  Entries to be used are randomly selected at the
// initialization time with equal probabilities to make the expected
// number of stored records equal to the specified parameter.  If the
// number of inputs is less than averageNumRecordsToUse, all input
// records are used.  Setting averageNumRecordsToUse<=0 disables the
// feature (the default).
//
// The fire() method returns one of the stored entries.  All the stored
// entries are sampled with equal probabilities.
//
// See StoppedParticleReactionGun_module.cc for an example of use.
//
// Andrei Gaponenko, 2014

#ifndef RootTreeSampler_hh
#define RootTreeSampler_hh

#include <vector>

#include "fhiclcpp/ParameterSet.h"

#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandFlat.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"

namespace mu2e {

  template<class Record>
  class RootTreeSampler {
  public:

    RootTreeSampler(art::RandomNumberGenerator::base_engine_t& engine,
                    const fhicl::ParameterSet& pset,
                    // a rudimentary protection against wrong input trees, which
                    // can cause memory overwrites
                    unsigned nBranchLeaves);

    const Record& fire() { return records_.at(randFlat_.fireInt(records_.size())); }

    typename std::vector<Record>::size_type
    numRecords() const { return records_.size(); }

  private:
    CLHEP::RandFlat randFlat_;
    std::vector<Record> records_;

    typedef std::vector<std::string> Strings;

    long countInputRecords(const art::ServiceHandle<art::TFileService>& tfs,
                           const Strings& files,
                           const std::string& treeName);
  };

}

// The user interface ends here.
//================================================================
// The implementation follows.

#include "cetlib/exception.h"

#include "TTree.h"
#include "TFile.h"

#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"

namespace mu2e {

  template<class Record> RootTreeSampler<Record>::
  RootTreeSampler(art::RandomNumberGenerator::base_engine_t& engine,
                  const fhicl::ParameterSet& pset,
                  unsigned nBranchLeaves)
    : randFlat_(engine)
  {
    const auto inputFiles(pset.get<std::vector<std::string> >("inputFiles"));
    const auto treeName(pset.get<std::string>("treeName"));
    const auto branchName(pset.get<std::string>("branchName"));
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
          std::cout<<"RootTreeSampler: recordUseFraction = "<<1
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

      const Long64_t nTreeEntries = nt->GetEntries();
      if(verbosityLevel > 0) {
        std::cout<<"RootTreeSampler: reading "<<nTreeEntries
                 <<" entries.  Tree "<<treeName
                 <<", file "<<infile->GetName()
                 <<std::endl;
      }

      Record rr;
      TBranch *bb = nt->GetBranch(branchName.c_str());
      if(!bb) {
        throw cet::exception("BADINPUT")<<"RootTreeSampler: Could not get branch \""<<branchName
                                        <<"\" in tree \""<<treeName
                                        <<"\" from file \""<<infile->GetName()
                                        <<"\"\n";
      }

      if(unsigned(bb->GetNleaves()) != nBranchLeaves) {
        throw cet::exception("BADINPUT")<<"RootTreeSampler: wrong number of leaves: expect "
                                        <<nBranchLeaves<<", but branch \""<<branchName
                                        <<"\", tree \""<<treeName
                                        <<"\" in file \""<<infile->GetName()
                                        <<"\" has "<<bb->GetNleaves()
                                        <<"\n";
      }

      bb->SetAddress(&rr);

      records_.reserve(records_.size()
                       // Add "mean + 3sigma": do not re-allocate in most cases.
                       + recordUseFraction*nTreeEntries
                       + 3*recordUseFraction*sqrt(double(nTreeEntries)));

      for(Long64_t i=0; i<nTreeEntries; ++i) {
        bb->GetEntry(i);

        if((averageNumRecordsToUse <= 0) || (randFlat_.fire() < recordUseFraction)) {
          records_.emplace_back(rr);
        }
      }
    }
  }

  //================================================================
  template<class Record>
  long RootTreeSampler<Record>::countInputRecords(const art::ServiceHandle<art::TFileService>& tfs,
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
