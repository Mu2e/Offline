#ifndef BTrkHelper_FileFinder_hh
#define BTrkHelper_FileFinder_hh

#include "BTrk/BaBar/FileFinderInterface.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"

#include "fhiclcpp/ParameterSet.h"

#include <string>

namespace mu2e {

  class FileFinder : public FileFinderInterface {

  public:
    FileFinder( fhicl::ParameterSet const& pset);

    std::string matElmDictionaryFileName() const override;
    std::string matIsoDictionaryFileName() const override;
    std::string matMtrDictionaryFileName() const override;

    std::string findFile( std::string const& ) const override;

  private:

    mutable ConfigFileLookupPolicy policy_;

    std::string elementsBaseName_;
    std::string isotopesBaseName_;
    std::string materialsBaseName_;

  };

}

#endif /* btrkHelper_FileFinder_hh */
