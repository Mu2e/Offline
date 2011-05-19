#ifndef LTrackerGeom_LayerInfo_hh
#define LTrackerGeom_LayerInfo_hh
//
// Information about a Layer.  Used by LTrackerMaker
// to construct an LTracker.
//
//
// $Id: LayerInfo.hh,v 1.5 2011/05/19 21:53:36 wb Exp $
// $Author: wb $
// $Date: 2011/05/19 21:53:36 $
//
// Original author Rob Kutschke
//

namespace mu2e {

  class LayerInfo{

  public:

    // The straws with adjacent cathode pads are conductive.
    // The straws without adjacent cathode pads are non-conductive.
    enum Stype {conductive, nonconductive, undefined};

    LayerInfo():
      _nStraws(-1),
      _strawType(undefined){
    }

    LayerInfo( int nStraws,
               Stype strawType
               ):
      _nStraws(nStraws),
      _strawType(strawType){
    }

    // Use compiler-generated copy c'tor, copy assignment, and d'tor

  int nStraws() const { return _nStraws; }
  Stype strawType() const { return _strawType; }

  private:

    int _nStraws;
    Stype _strawType;

  };

}  //namespace mu2e

#endif /* LTrackerGeom_LayerInfo_hh */
