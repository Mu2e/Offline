#ifndef DisplayUtilsConfig_hh
#define DisplayUtilsConfig_hh

#include  "Offline/TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eMainWindow.h"

#include "fhiclcpp/types/Atom.h"

namespace mu2e {

  namespace DisplayUtilsConfig {

    // Configuration of the "show" table within the top level configuration.
    struct ShowConfig {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<bool> showbuilding{Name("showBuilding"), Comment("set false to remove building"),false};
      fhicl::Atom<bool> showCRV{Name("showCRV"), Comment("set false if you just want to see DS"),false};
      fhicl::Atom<bool> showDSOnly{Name("showDSOnly"), Comment(""),true};
      fhicl::Atom<bool> showInsidePS{Name("showInsidePS"), Comment(""),false};
    };

    // Function to populate the Show struct from the validated config.
    inline TEveMu2eMainWindow::GeomOptions convert( ShowConfig const& c ){
    return TEveMu2eMainWindow::GeomOptions(c.showbuilding(),
                       c.showCRV(),
                       c.showDSOnly(),
                       c.showInsidePS() );
    }

  }

}

#endif 
