#include "BtrkHelper/inc/FileFinder.hh"

mu2e::FileFinder::FileFinder( fhicl::ParameterSet const& pset ):
  policy_(),
  elementsBaseName_(pset.get<std::string>("elements")),
  isotopesBaseName_(pset.get<std::string>("isotopes")),
  materialsBaseName_(pset.get<std::string>("materials")){
}

std::string mu2e::FileFinder::matElmDictionaryFileName() const{
  return findFile(elementsBaseName_);
}

std::string mu2e::FileFinder::matIsoDictionaryFileName() const{
  return findFile(isotopesBaseName_);
}

std::string mu2e::FileFinder::matMtrDictionaryFileName() const{
  return findFile(materialsBaseName_);
}

std::string mu2e::FileFinder::findFile( std::string const& path ) const{
  return policy_( path );
}
