#ifndef ROOTtools_artProductSizes_RootSizeOnDisk_hh
#define ROOTtools_artProductSizes_RootSizeOnDisk_hh
//
// Collect information about the disk space used by the top tier and second
// tier objects inside an art format root event-data file.
//
// This is known to work for root 5.34.*
//
//
// The available information is:
//
//    filename  - The name of the root file
//    size      - The size on disk according to the file system; equivalent to ls -l
//    sum       - The sum of the sizes of all top tier TTree objects.
//                Some day this may include the sizes of top tier objects that are not TTrees.
//    fraction  - The ratio of sum()/size()
//    contents  - Detailed information about each TTree, including the size of each branch.
//

#include <ostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"

namespace mu2e {

  class RootSizeOnDisk {

  public:
    RootSizeOnDisk ( std::string const& aFileName, TFile* aFile );

    // Information about a top tier or second tier object in an art format root event-data file.
    class Record {

    public:

      Record ( std::string const& aname,
               std::string const& aclassName,
               Long64_t           asize=0,
               Long64_t           acount=1,
               double             afraction=0. );
      bool operator< ( Record const& rhs) const{
        return (name() < rhs.name());
      }

      bool isTree() const{
        return (className_ == "TTree");
      }

      bool isTKey() const{
        return (className_ == "TKey");
      }

      std::string const&         name()      const { return name_; }
      std::string const&         className() const { return className_; }
      Long64_t                   size()      const { return size_; }
      Long64_t                   count()     const { return count_; }
      double                     fraction()  const { return fraction_; }
      std::vector<Record> const& contents()  const { return contents_; }

      double sizePerEvent() const;

      // Modifiers
      void size ( Long64_t s ) { size_     = s;}
      void count( Long64_t c ) { count_    = c;}
      void fraction( double f) { fraction_ = f;}
      void contents( std::vector<Record>& c ) { contents_.swap(c); }

    private:
      std::string name_;       // Name of the object
      std::string className_;  // Classname, either TTree or TKey
      Long64_t    size_;       // Total space on disk occupied by this object
      Long64_t    count_;      // If a TTree, the number of events in the TTree; if a TKey then it is 1
      double      fraction_;   // Fraction of the space occupied by objects at this level that is occupied by this object

      std::vector<Record> contents_;

    }; // end RootSizeOnDisk::Record

    //typedef std::vector<RootSizeOnDisk::Record> Records_t;
    typedef std::vector<Record> Records_t;

    std::string const& filename()  const { return fileName_; }
    Long64_t           size()      const { return size_;     }
    Long64_t           sum()       const { return sum_;      }
    double             fraction()  const { return fraction_; }
    Records_t const&   contents()  const { return contents_; }

    void print ( std::ostream& os, double minimumFraction ) const;

  private:

    std::string fileName_;
    Long64_t    size_;
    Long64_t    sum_;
    double      fraction_;
    Records_t   contents_;

    void fillLevel2( Record&, TTree* );

  };

  // Compare two RootSizeOnDisk::Record objects to select the one with the larger size.
  bool greaterBySize( RootSizeOnDisk::Record const& lhs, RootSizeOnDisk::Record const& rhs );

} // namespace mu2e

#endif /* ROOTtools_artProductSizes_RootSizeOnDisk_hh */
