#ifndef Analyses_GeneratorSummaryHistograms_hh
#define Analyses_GeneratorSummaryHistograms_hh
//
// Make histograms summarizing the information in the event generator.
//
// $Id: GeneratorSummaryHistograms.hh,v 1.2 2013/09/08 01:30:05 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/09/08 01:30:05 $
//
// Contact person Rob Kutschke
//

#include "MCDataProducts/inc/GenParticleCollection.hh"

#include "art_root_io/TFileDirectory.h"

#include <string>

// Forward declarations

class TH1F;
class TH2F;

namespace mu2e {

  class DetectorSystem;

  class GeneratorSummaryHistograms{

  public:

    GeneratorSummaryHistograms();
    // Accept compiler generated d'tor.  Class is not copyable; see private section.

    // Book histograms at the root TFileDirectory for the current module.
    void book( );

    // Book histograms in the subdirectory, given by the relativePath; that path is
    // relative to the root TFileDirectory for the current module.
    void book( std::string const& relativePath );

    // Book histograms in the specified TFileDirectory.
    void book( art::TFileDirectory& tfdir );

    void fill( GenParticleCollection const& sims );

  private:

    // Not copyable or assignable.  These will not be implemented.
    GeneratorSummaryHistograms ( GeneratorSummaryHistograms const& rhs );
    GeneratorSummaryHistograms& operator=(GeneratorSummaryHistograms const& rhs);

    DetectorSystem const* detSys_;

    TH1F* hMultiplicity_;
    TH1F* hgenId_;
    TH1F* hp_;
    TH1F* hpt_;
    TH1F* hke_;
    TH1F* hcz_;
    TH1F* hphi_;
    TH1F* hradius_;
    TH1F* hz_;
    TH1F* htime_;
    TH2F* hxy_;
    TH2F* hrz_;

  }; // end class Diagnostics G4

} // end namespace mu2e

#endif /* Analyses_GeneratorSummaryHistograms_hh */
