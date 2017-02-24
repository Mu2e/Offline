//
// Scan the Geant4 logical volume store and look for
// logical volume names that are duplicated.
//

// Mu2e includes
#include "Mu2eG4/inc/DuplicateLogicalVolumeChecker.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"

#include "cetlib_except/exception.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <ostream>

mu2e::DuplicateLogicalVolumeChecker::DuplicateLogicalVolumeChecker():
  duplicates_(){

  G4LogicalVolumeStore*  vstore = G4LogicalVolumeStore::GetInstance();
  if ( vstore == 0 ){
    throw cet::exception("GEOM")
      << __func__
      << "  could not access G4LogicalVolumeStore pointer.\n";
  }
  storeSize_=vstore->size();

  // Temporary to accumulate the information.
  std::map<std::string,int> lnames;
  for ( auto i : *vstore ){
    ++lnames[i->GetName()];
  }

  // No duplicates.
  if ( vstore->size() == lnames.size() ){
    return;
  }

  // Capacity is positive by construction
  duplicates_.reserve( vstore->size()-lnames.size());

  // Prepare the summary information from the temporary info.
  for ( auto const& i : lnames ){
    if ( i.second > 1 ){
      duplicates_.push_back( Info(i.first,i.second) );
    }
  }

}

void mu2e::DuplicateLogicalVolumeChecker::print( std::ostream& out ) const{

  if ( duplicates_.empty() ){
    out << "The G4LogicalVolumeStore has no duplicate names." << std::endl;
    out << "Size of the G4LogicalVolumeStore: " << storeSize_ << std::endl;
    return;
  }

  int widthName{0};
  int maxN{1};
  for ( auto const& i : duplicates_ ){
    widthName = std::max( widthName, int(i.name.size()) );
    maxN      = std::max( maxN,      i.n );
  }

  widthName += 3;
  int widthN = static_cast<int>( log10(double(maxN)))+3;
  widthN = std::min( widthN, 7);

  std::string gap("   ");

  // Print the header
  out << "\nThe G4LogicalVolumeStore has duplicate names: " << std::endl;
  out << std::setw(widthN) << "Count"
      << gap
      << "Logical Volume Name"
      << std::endl;
  for ( auto const& i : duplicates_ ){
    out << std::setw(widthN)    << i.n
        << gap << i.name
        << std::endl;
  }

}

bool mu2e::DuplicateLogicalVolumeChecker::hasForbiddenNames( std::ostream& out,
                                                             std::vector<std::string>const& forbidden,
                                                             bool throwOnError ) const{

  int nbad{0};
  for ( auto const& info : duplicates_ ){
    auto i = find( forbidden.begin(), forbidden.end(), info.name);
    if ( i != forbidden.end() ) {
      out << "Logical volume name on forbidden list is duplicate: " << *i << std::endl;
      ++nbad;
    }
  }

  if ( nbad != 0 ){
    if  ( throwOnError ) {
      throw cet::exception("GEOM")
        << __func__
        << " there are duplicate names in G4LogicalVolumeStore.\n"
        << " Some of these are on the forbiddent list.\n";
    }
  }

  return (nbad!=0);
}

mu2e::DuplicateLogicalVolumeChecker::Info::Info ( std::string const& aname, int an):
  name(aname), n(an){
}
