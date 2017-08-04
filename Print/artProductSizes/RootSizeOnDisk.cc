//
// Collect information about the disk space used by the top tier and second
// tier objects inside an art format root event-data file.
//

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <set>

#include "Print/artProductSizes/RootSizeOnDisk.hh"
#include "Print/artProductSizes/rootFileSizeTools.hh"

#include "TFile.h"
#include "TList.h"
#include "TKey.h"
#include "TTree.h"
#include "TBranch.h"
#include "TObjArray.h"

#include "boost/filesystem.hpp"
#include "boost/format.hpp"

using namespace std;

namespace {

  // Helper struct used to pull information out of a TList.
  struct MyObjects {

    MyObjects( set<mu2e::RootSizeOnDisk::Record>& akeys):
      keys(akeys){
    }

    bool operator()(TObject *aObj) {
      TKey* key = (TKey*) aObj;
      keys.insert( mu2e::RootSizeOnDisk::Record( key->GetName(), key->GetClassName()) );
      return true;
    }

    set<mu2e::RootSizeOnDisk::Record>& keys;
  };

} // end anonymous namespace

mu2e::RootSizeOnDisk::Record::Record ( std::string const& aname,
                                       std::string const& aclassName,
                                       Long64_t           asize,
                                       Long64_t           acount,
                                       double             afraction ):
  name_(aname),
  className_(aclassName),
  size_(asize),
  count_(acount),
  fraction_(afraction){
}

double mu2e::RootSizeOnDisk::Record::sizePerEvent() const{
  return ( count_ > 0 ) ? size_/double(count_) : size_;
}

bool mu2e::greaterBySize( mu2e::RootSizeOnDisk::Record const& lhs, mu2e::RootSizeOnDisk::Record const& rhs ){
  return ( lhs.size() > rhs.size() );
}

void mu2e::RootSizeOnDisk::print( std::ostream& os, double minimumFraction ) const{

  os << "\nSize on disk for the file: " << filename() << "\n"
     << "Size on disk (bytes):      " << size()
     << endl;
  os << setw(13) << "Size"
     << setw(13) << "Size/Entry"
     << setw(14) << "Entries"
     << setw(12) << "Fraction"
     << "  TTree/TKey Name"

     << endl;
  os << setw(13) << "(bytes)"
     << setw(13) << "(bytes)"
     << endl;

  for ( RootSizeOnDisk::Record const& key : contents() ){
    if ( key.isTree() || key.isTKey() ){
      os << setw(13) << key.size();
      if ( key.isTree() ) {
        os << boost::format("%13.1f") % key.sizePerEvent()
           << boost::format("%14d") % key.count();
      } else{
        os << setw(13) << "--"
           << setw(14) << "--";
      }
      os << boost::format("%12.3f") % key.fraction() << "  "
         << key.name()
         << endl;
    } else {
      os << setw(13) << key.size() << " "
         << boost::format("%10.3f") % key.fraction() << " "
         << key.name()
         << "  (skipped because not a TTree or a TKey; it is a"
         << key.className() << ")"
         << endl;
    }
  }
  os << "----------------------------------------------------\n"
     << setw(13) << sum()
     << setw(29) << " "
     << boost::format("%10.3f") % fraction() << "  "
     << "Total\n"
     << endl;

  os << "Details for each TTree that occupies more than the fraction "
     << minimumFraction
     << " of the size on disk.\n"
     << endl;

  for ( RootSizeOnDisk::Record const& key : contents() ){
    if ( key.isTree() && (key.fraction() > minimumFraction ) ){
      os << "\nDetails for branch: "  << key.name() << endl;
      os << "Number of entries:  " << key.count() << endl;
      os << setw(12) << "Size"
         << setw(15) << "Size/Entry"
         << setw(10) << "   Fraction"
         << "  Data Product Name"
         << endl;
      os << setw(12) << "(bytes)"
         << setw(15) << "(bytes)\n"
         << endl;

      Long64_t sum(0);
      for ( auto const& branch : key.contents() ){
        sum += branch.size();
        os << setw(12) << branch.size() << " "
           << boost::format("%14.1f") % branch.sizePerEvent()  << " "
           << boost::format("%10.3f") % branch.fraction()      << "  "
           << branch.name()
           << endl;
      }
      double ratio = double(sum)/double(key.size());
      double sizePerEvent = (key.count()>0) ? double(sum)/double(key.count()) : sum;
      os << "--------------------------------------\n"
         << setw(12) << sum << " "
         << boost::format("%14.1f") % sizePerEvent  << " "
         << boost::format("%10.3f") % ratio
         << "  Total\n"
         << endl;
    }
  }
}

mu2e::RootSizeOnDisk::RootSizeOnDisk ( std::string const& aFileName, TFile* file ):
  fileName_(aFileName),
  size_(0),
  sum_(0),
  fraction_(0){

  // File size on disk, in bytes.
  size_ = boost::filesystem::file_size( fileName_.c_str() );

  // Extract info about top level objects.
  // There are usually Multiple cycles of these objects; we only want each name once.
  set<RootSizeOnDisk::Record> topKeys;
  TList* keys = file->GetListOfKeys();
  TIter iter(keys);
  for_each( iter.Begin(), TIter::End(), MyObjects(topKeys) );

  // Copy to a vector since we will want to sort them by size.
  contents_.assign( topKeys.begin(), topKeys.end() );

  // Compute sizes of each top level TTree and TKey.
  for ( vector<RootSizeOnDisk::Record>::iterator i=contents_.begin(), e=contents_.end(); i!=e; ++i){
    RootSizeOnDisk::Record & key(*i);
    if ( key.isTree() ){
      TTree * tree;
      file->GetObject( key.name().c_str(), tree);
      Long64_t size = sizeOnDisk(tree);
      sum_ += size;
      double f = double(size)/double(size_);
      key.size(size);
      key.count(tree->GetEntries());
      key.fraction(f);
      fillLevel2( key, tree );
    } else if ( key.isTKey() ){
      TKey * tkey = file->FindKey( key.name().c_str());
      Long64_t size = tkey->GetNbytes();
      sum_ += size;
      double f = double(size)/double(size_);
      key.size(size);
      key.fraction(f);
    }
  }

  // Sort by decreasing size.
  sort ( contents_.begin(), contents_.end(), greaterBySize );

  fraction_ = double(sum_)/double(size_);

}

// For each TTree Record, fill the information about the branches.
void mu2e::RootSizeOnDisk::fillLevel2( RootSizeOnDisk::Record& key, TTree* tree){

  TObjArray* branches = tree->GetListOfBranches();
  size_t n = branches->GetEntries();

  Records_t branchInfo;

  for( size_t i = 0; i < n; ++i ) {
    TBranch *subbr = dynamic_cast<TBranch*>(branches->At(i));
    Long64_t size =  sizeOnDisk(subbr,true);
    double f = double(size)/double(key.size());
    branchInfo.emplace_back( subbr->GetName(), "TBranch", size, key.count(), f );
  }

  sort( branchInfo.begin(), branchInfo.end(), greaterBySize);
  key.contents( branchInfo);

}
