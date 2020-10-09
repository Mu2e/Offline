#ifndef ROOTtools_artProductSizes_rootFileSizeTools_hh
#define ROOTtools_artProductSizes_rootFileSizeTools_hh
//
// Utilties for finding the size on disk of TTrees and TBranches
// Adapted from code given to us by Philippe Canal, pcanal@fnal.gov.
//
// This is known to work for root 5.34.*
//
//

#include "TBranch.h"
#include "TObjArray.h"
#include "TTree.h"

namespace mu2e {

  Long64_t GetBasketSize(TObjArray * branches, bool ondisk, bool inclusive);
  Long64_t GetBasketSize(TBranch * b, bool ondisk, bool inclusive);
  Long64_t GetTotalSize( TBranch * br, bool ondisk, bool inclusive );
  Long64_t GetTotalSize( TObjArray * branches, bool ondisk );
  Long64_t GetTotalSize(TTree *t, bool ondisk);
  Long64_t sizeOnDisk(TTree *t);
  Long64_t sizeOnDisk(TBranch *branch, bool inclusive);
  void     printBranchSummary(TBranch *br);
  void     printTreeSummary(TTree *t);

} // namespace mu2e

#endif /* ROOTtools_artProductSizes_rootFileSizeTools_hh */
