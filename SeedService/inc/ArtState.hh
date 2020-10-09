#ifndef SeedService_ArtState_hh
#define SeedService_ArtState_hh
//
// Describe the current state of art processing, as understood by the SeedService.
//
//
// Contact person Rob Kutschke
//

#include <string>

namespace mu2e {

  namespace SeedServiceHelper {

    struct ArtState {

      enum state_type { unDefined, inConstructor, inBeginRun };

      ArtState():
        state(unDefined),
        currentModuleLabel(){
      }

      // Accept compiler written d'tor, copy c'tor and copy assignment.

      void set( state_type astate, std::string const& label ){
        state = astate;
        currentModuleLabel = label;
      }

      void clear(){
        state = unDefined;
        currentModuleLabel = "";
      }

      state_type state;
      std::string currentModuleLabel;

    }; // end ArtState

  } // end namespace SeedServiceHelper

} // end namespace mu2e

#endif /* SeedService_ArtState_hh */
