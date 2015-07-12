#ifndef Mu2eG4_DuplicateLogicalVolumeChecker_hh
#define Mu2eG4_DuplicateLogicalVolumeChecker_hh

//
// Scan the Geant4 logical volume store and look for
// logical volume names that are duplicated.
//

#include <iosfwd>
#include <vector>
#include <string>

namespace mu2e {


  class DuplicateLogicalVolumeChecker {
  public:

    DuplicateLogicalVolumeChecker();

    // Name and count of logical volume names that appear more than once.
    struct Info {
      Info ( std::string const& aname, int an);
      std::string name;
      int         n;
    };

    // If true, there are no duplicated logical volume names.
    bool ok()            const { return duplicates_.empty(); }

    // If true, there are duplicated logical volume names.
    bool hasDuplicates() const { return !duplicates_.empty(); }

    // If true, the set of duplicated names includes one or more
    // of the forbidden names, supplied as an argument.
    bool hasForbiddenNames( std::ostream& out,
                            std::vector<std::string>const& fobidden,
                            bool throwIfError = true ) const;

    // Full information about the duplicates.
    auto const& duplicates() const { return duplicates_; }

    // Print the full info.
    void print( std::ostream& out) const;

  private:

    // Size of the logical volume store.
    int            storeSize_ = -1;

    // Information about duplicate names.
    std::vector<Info> duplicates_;

  };

} // end namespace mu2e

#endif /* Mu2eG4_DuplicateLogicalVolumeChecker_hh */
