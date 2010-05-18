#ifndef STRAWPLACER_HH
#define STRAWPLACER_HH
//
// Class to place one straw within the tracker mother volume.
//
// $Id: StrawPlacer.hh,v 1.2 2010/05/18 20:28:41 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/05/18 20:28:41 $
//
// Original author Rob Kutschke
//
// This implementation does not know about segmentation of the 
// tracker mother volume. Segmenation should be added to improve 
// navigation speed.
// 

#include <string>

class G4LogicalVolume;

namespace mu2e {

  class Straw;

  class StrawPlacer{

  public:

    StrawPlacer( std::string basename
                 , G4LogicalVolume* logical
                 , G4LogicalVolume* motherLogical
                 );
    ~StrawPlacer();

    void operator() ( const Straw& s );

  private:
    std::string      _basename;
    G4LogicalVolume* _logical;
    G4LogicalVolume* _motherLogical;

  };

}

#endif
