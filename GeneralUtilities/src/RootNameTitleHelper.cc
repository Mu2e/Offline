//
// Helper class to construct a series of similar
// names and titles for root objects.  See the header
// for details.
//
//

#include <iostream>
#include <iomanip>
#include <sstream>

#include "Offline/GeneralUtilities/inc/RootNameTitleHelper.hh"

using namespace std;

RootNameTitleHelper::RootNameTitleHelper(
                                         std::string const& name_base,
                                         std::string const& title_base,
                                         int id,
                                         int pad ){

  // Form the name.
  ostringstream name;
  name.fill('0');
  if ( pad > 0 ){
    name << name_base << setw(pad) << id;
  } else{
    name << name_base << id;
  }
  _name = name.str();

  // Form the title.
  ostringstream title;
  title << title_base << id;
  _title = title.str();


}
