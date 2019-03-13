#ifndef BTrkHelper_FileFinder_hh
#define BTrkHelper_FileFinder_hh

#include "BTrk/BaBar/FileFinderInterface.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"

#include "fhiclcpp/ParameterSet.h"

#include <string>

namespace mu2e {

  class FileFinder : public FileFinderInterface {

  public:
    FileFinder(std::string elementsBaseName,
	       std::string isotopesBaseName,
	       std::string materialsBaseName):
      elementsBaseName_(elementsBaseName),
      isotopesBaseName_(isotopesBaseName),
      materialsBaseName_(materialsBaseName){}

    std::string matElmDictionaryFileName() const override {
      return findFile(elementsBaseName_);
    }

    std::string matIsoDictionaryFileName() const override {
      return findFile(isotopesBaseName_);
    }

    std::string matMtrDictionaryFileName() const override {
      return findFile(materialsBaseName_);
    }

    std::string findFile( std::string const& path ) const override {
      return policy_( path );
    }

  private:

    // mutable because BTrk want to hold this with
    // a const pointer, but this needs to alter its internal state
    mutable ConfigFileLookupPolicy policy_;

    std::string elementsBaseName_;
    std::string isotopesBaseName_;
    std::string materialsBaseName_;

  };

}

#endif /* btrkHelper_FileFinder_hh */
