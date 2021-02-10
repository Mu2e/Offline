#ifndef Mu2eKinKal_KKFileFinder_hh
#define Mu2eKinKal_KKFileFinder_hh

#include "KinKal/MatEnv/FileFinderInterface.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include <string>

namespace mu2e {

  class KKFileFinder : public MatEnv::FileFinderInterface {

  public:
    KKFileFinder(std::string elementsBaseName,
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

#endif /* btrkHelper_KKFileFinder_hh */
