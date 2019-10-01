//
// Class which extracts informayion from the framework event objects to build the event display shapes (e.g. tracks, straws, support structures).
//
// $Id: DataInterface.h,v 1.33 2014/02/22 01:52:18 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2014/02/22 01:52:18 $
//
// Original author Ralf Ehrlich
//

#ifndef EventDisplay_src_DataInterface_h
#define EventDisplay_src_DataInterface_h

#include "CLHEP/Vector/ThreeVector.h"
#include "EventDisplay/src/ContentSelector.h"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "Mu2eBTrk/inc/ParticleInfo.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "art/Framework/Principal/Event.h"
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
class Cylinder;
class Cone;
class ComponentInfo;
class ContentSelector;
class EventDisplayFrame;

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
  EventDisplayFrame  *_mainframe;  //points to the EventDisplayFrame object
                                   //ROOT needs a bare pointer to this object when dealing
                                   //with context menus (the function which gets called
                                   //from the context menu belongs to this object)
  std::list<boost::shared_ptr<VirtualShape> >       _components;
  std::map<int, boost::shared_ptr<Straw> >          _straws;
  std::map<int, boost::shared_ptr<VirtualShape> >   _crystals;
  std::map<int, boost::shared_ptr<Cube> >           _crvscintillatorbars;
  std::vector<boost::shared_ptr<Straw> >        _hits;
  std::vector<boost::shared_ptr<VirtualShape> > _crystalhits;
  std::vector<boost::shared_ptr<Cylinder> >     _driftradii;
  std::vector<boost::shared_ptr<Track> >        _tracks;
  std::vector<boost::shared_ptr<Cube> >         _crvhits;
  std::vector<boost::shared_ptr<VirtualShape> > _supportstructures;
  std::vector<boost::shared_ptr<VirtualShape> > _otherstructures;
  std::vector<boost::shared_ptr<VirtualShape> > _mbsstructures;
  std::vector<boost::shared_ptr<Cone> > _mecostylepastructures;
  CLHEP::Hep3Vector _detSysOrigin;
  timeminmax        _hitsTimeMinmax, _tracksTimeMinmax;
  spaceminmax       _trackerMinmax, _targetMinmax, _calorimeterMinmax, _tracksMinmax;
  int               _numberHits, _numberCrystalHits;
  bool              _showUnhitStraws, _showUnhitCrystals;
  unsigned int _minPoints;
  double _minTime;
  double _maxTime;
  double _minMomentum;
  bool _showElectrons;
  bool _showMuons;
  bool _showGammas;
  bool _showNeutrinos;
  bool _showNeutrons;
  bool _showOthers;
  mu2e::StrawHitFlag _hitFlagSetting;

  std::unique_ptr<mu2e::ParticleInfo> _particleInfo;

  void createGeometryManager();
  void removeAllComponents();
  void removeNonGeometryComponents();
  void findBoundaryT(timeminmax &m, double t);
  void findBoundaryP(spaceminmax &m, double x, double y, double z);
  void resetBoundaryT(timeminmax &m);
  void resetBoundaryP(spaceminmax &m);
  void toForeground();
  void findTrajectory(boost::shared_ptr<ContentSelector> const &contentSelector,
                      boost::shared_ptr<Track> const &track, const cet::map_vector_key &id,
                      double timeOffset,
                      const ContentSelector::trackInfoStruct &trackInfo);
  struct trajectoryStruct
  {
    CLHEP::Hep3Vector v;
    double t;
    trajectoryStruct() {t=NAN;}
  };

  public:
  DataInterface(EventDisplayFrame *mainframe);
  virtual ~DataInterface();

  void startComponents();
  void updateComponents(double time, boost::shared_ptr<ContentSelector> contentSelector);
  void fillGeometry();
  void fillEvent(boost::shared_ptr<ContentSelector> const &contentSelector, const mu2e::SimParticleTimeOffset &timeOffsets);
  void makeSupportStructuresVisible(bool visible);
  void makeOtherStructuresVisible(bool visible);
  void makeCrvScintillatorBarsVisible(bool visible);
  void makeStrawsVisibleBeforeStart(bool visible);
  void makeCrystalsVisibleBeforeStart(bool visible);
  void makeMuonBeamStopStructuresVisible(bool visible);
  void makeMecoStyleProtonAbsorberVisible(bool visible);
  void useHitColors(bool hitcolors, bool whitebackground);
  void useTrackColors(boost::shared_ptr<ContentSelector> const &contentSelector, bool trackcolors, bool whitebackground);
  void getFilterValues(unsigned int &minPoints, double &minTime, double &maxTime, double &minMomentum,
                       bool &showElectrons, bool &showMuons, bool &showGammas, 
                       bool &showNeutrinos, bool &showNeutrons, bool &showOthers,
                       mu2e::StrawHitFlag &hitFlagSetting);
  void setFilterValues(unsigned int minPoints, double minTime, double maxTime, double minMomentum,
                       bool showElectrons, bool showMuons, bool showGammas, 
                       bool showNeutrinos, bool showNeutrons, bool showOthers,
                       mu2e::StrawHitFlag hitFlagSetting);
  int getNumberHits() {return _numberHits;}
  int getNumberCrystalHits() {return _numberCrystalHits;}

  timeminmax getHitsTimeBoundary(); 
  timeminmax getTracksTimeBoundary();

  spaceminmax getTrackerBoundary() {return _trackerMinmax;}
  spaceminmax getTargetBoundary() {return _targetMinmax;}
  spaceminmax getCalorimeterBoundary() {return _calorimeterMinmax;}
  spaceminmax getTracksBoundary() {return _tracksMinmax;}
  spaceminmax getSpaceBoundary(bool useTarget, bool useCalorimeter, bool useTracks);
};

}
#endif /* EventDisplay_src_DataInterface_h */
