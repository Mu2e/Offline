//
// Class which extracts informayion from the framework event objects to build the event display shapes (e.g. tracks, straws, support structures).
//
// $Id: DataInterface.h,v 1.17 2011/08/07 18:21:06 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2011/08/07 18:21:06 $
//
// Original author Ralf Ehrlich
//

#ifndef EventDisplay_src_DataInterface_h
#define EventDisplay_src_DataInterface_h

#include "CLHEP/Vector/ThreeVector.h"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "art/Framework/Core/Event.h"
#include "boost/shared_ptr.hpp"
#include <TObject.h>
#include <list>
#include <map>

class TGeoManager;
class TGeoVolume;
class TGMainFrame;

namespace mu2e_eventdisplay
{

class VirtualShape;
class Track;
class Straw;
class Cube;
class Sphere;
class ComponentInfo;
class ContentSelector;

class DataInterface
{
  DataInterface();
  DataInterface(const DataInterface &);
  DataInterface& operator=(const DataInterface &);

  public:
  struct timeminmax
  {
    double mint, maxt;
  };
  struct spaceminmax
  {
    double minx, maxx;
    double miny, maxy;
    double minz, maxz;
  };

  private:
  art::Event    *_event;
  TGeoManager   *_geometrymanager; //bare pointer needed since ROOT manages this object
  TGeoVolume    *_topvolume;       //bare pointer needed since ROOT manages this object
  const TObject *_mainframe;       //points to the EventDisplayFrame object
                                   //ROOT needs a bare pointer to this object when dealing
                                   //with context menus (the function which gets called
                                   //from the context menu belongs to this object)
  std::list<boost::shared_ptr<VirtualShape> >   _components;
  std::map<int, boost::shared_ptr<Straw> >      _straws;
  std::map<int, boost::shared_ptr<Cube> >       _crystals;
  std::vector<boost::shared_ptr<Straw> >        _hits;
  std::vector<boost::shared_ptr<Cube> >         _crystalhits;
  std::vector<boost::shared_ptr<Sphere> >       _spheres;
  std::vector<boost::shared_ptr<Track> >        _tracks;
  std::vector<boost::shared_ptr<VirtualShape> > _supportstructures;
  std::vector<boost::shared_ptr<VirtualShape> > _otherstructures;
  double            _xOffset, _zOffset, _zOffsetDS;
  CLHEP::Hep3Vector _mu2eOriginInWorld;
  timeminmax        _hitsTimeMinmax, _tracksTimeMinmax;
  spaceminmax       _trackerMinmax, _targetMinmax, _calorimeterMinmax, _tracksMinmax;
  int               _numberHits, _numberCrystalHits;
  bool              _showUnhitStraws, _showUnhitCrystals;

  void createGeometryManager();
  void removeNonGeometryComponents();
  void removeAllComponents();
  void findBoundaryT(timeminmax &m, double t);
  void findBoundaryP(spaceminmax &m, double x, double y, double z);
  void resetBoundaryT(timeminmax &m);
  void resetBoundaryP(spaceminmax &m);
  void toForeground();
  void findTrajectory(const ContentSelector *contentSelector,
                      boost::shared_ptr<Track> track, int id,
                      double t1, double t2,
                      const mu2e::SimParticleCollection *simParticles,
                      const std::vector<int> &daughterVect);
  struct trajectoryStruct
  {
    CLHEP::Hep3Vector v;
    double t;
    trajectoryStruct() {t=NAN;}
  };

  public:
  DataInterface(const TGMainFrame *mainframe);
  virtual ~DataInterface();

  void startComponents();
  void updateComponents(double time);
  void fillGeometry();
  void fillEvent(const ContentSelector *contentSelector);
  void makeSupportStructuresVisible(bool visible);
  void makeOtherStructuresVisible(bool visible);
  void makeStrawsVisibleBeforeStart(bool visible);
  void makeCrystalsVisibleBeforeStart(bool visible);
  void useHitColors(bool hitcolors, bool whitebackground);
  void useTrackColors(bool trackcolors, bool whitebackground);
  int getNumberHits() {return _numberHits;}
  int getNumberCrystalHits() {return _numberCrystalHits;}

  timeminmax getHitsTimeBoundary() {return _hitsTimeMinmax;}
  timeminmax getTracksTimeBoundary() {return _tracksTimeMinmax;}
  spaceminmax getTrackerBoundary() {return _trackerMinmax;}
  spaceminmax getTargetBoundary() {return _targetMinmax;}
  spaceminmax getCalorimeterBoundary() {return _calorimeterMinmax;}
  spaceminmax getTracksBoundary() {return _tracksMinmax;}
  spaceminmax getSpaceBoundary(bool useTarget, bool useCalorimeter, bool useTracks);
};

}
#endif /* EventDisplay_src_DataInterface_h */
