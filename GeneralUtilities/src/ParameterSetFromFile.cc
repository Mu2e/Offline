//
//
// Read a .fcl file and form a parameter set object.
//

#include "GeneralUtilities/inc/ParameterSetFromFile.hh"

#include "fhiclcpp/intermediate_table.h"
#include "fhiclcpp/parse.h"
#include "fhiclcpp/make_ParameterSet.h"

#include "cetlib/filepath_maker.h"

#include <iostream>
#include <iostream>

mu2e::ParameterSetFromFile::
ParameterSetFromFile( std::string const& fileName ):
  _fileName(fileName)
{

  cet::filepath_lookup policy("FHICL_FILE_PATH");
  fhicl::intermediate_table tbl = fhicl::parse_document(_fileName, policy);
  _pSet = fhicl::ParameterSet::make(tbl);

}


void
mu2e::ParameterSetFromFile::printKeys( std::ostream& out ) const
{

  std::vector<std::string> const& keys = _pSet.get_names();
  out << "\nParameter set read from file: " << _fileName << std::endl;
  out << "Number of keys: " << keys.size() << std::endl;
  for ( std::vector<std::string>::const_iterator i=keys.begin();
        i != keys.end(); ++i ){
    out << " " << i-keys.begin() << " " << *i << std::endl;
  }

}
