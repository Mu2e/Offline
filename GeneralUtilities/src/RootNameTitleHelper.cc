//
// Helper class to construct a series of similar
// names and titles for root objects.  See the header
// for details.
//
// $Id: RootNameTitleHelper.cc,v 1.2 2010/05/18 21:15:44 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/05/18 21:15:44 $
//

#include <iostream>
#include <iomanip>
#include <sstream>

#include "GeneralUtilities/inc/RootNameTitleHelper.hh"

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
