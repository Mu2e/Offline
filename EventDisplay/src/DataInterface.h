//
// Class which extracts informayion from the framework event objects to build the event display shapes (e.g. tracks, straws, support structures).
//
// $Id: DataInterface.h,v 1.2 2011/01/29 02:14:20 ehrlich Exp $
// $Author: ehrlich $ 
// $Date: 2011/01/29 02:14:20 $
//
// Original author Ralf Ehrlich
//

#ifndef DATAINTERFACE_H
#define DATAINTERFACE_H

#include <TObject.h>
#include <list>
#include <map>
#include "FWCore/Framework/interface/Event.h"
#include "boost/shared_ptr.hpp"

class TGeoManager;
class TGeoVolume;
class TGMainFrame;

namespace mu2e_eventdisplay
{

class VirtualShape;
class Track;
class Straw;
class ComponentInfo;

class DataInterface
{
  DataInterface();
  DataInterface(const DataInterface &);
  DataInterface& operator=(const DataInterface &);

  struct minmax
  {
    double minx, maxx;
    double miny, maxy;
    double minz, maxz;
    double mint, maxt;
  };

  TGeoManager   *_geometrymanager; //bare pointer needed since ROOT manages this object
  TGeoVolume    *_topvolume;       //bare pointer needed since ROOT manages this object
  const TObject *_mainframe;       //points to the EventDisplayFrame object
                                   //ROOT needs a bare pointer to this object when dealing
                                   //with context menus (the function which gets called
                                   //from the context menu belongs to this object) 
  std::list<boost::shared_ptr<VirtualShape> >   _components;
  std::map<int, boost::shared_ptr<Straw> >      _straws;
  std::vector<boost::shared_ptr<Straw> >        _hits;
  std::vector<boost::shared_ptr<Track> >        _tracks;
  std::vector<boost::shared_ptr<VirtualShape> > _supportstructures;
  double        _z0;
  minmax        _hitsMinmax, _tracksMinmax;
  int           _numberHits;
  bool          _showUnhitStraws;

  void createGeometryManager();
  void removeNonGeometryComponents();
  void removeAllComponents();

  public:
  DataInterface(const TGMainFrame *mainframe);
  virtual ~DataInterface(); 

  const std::list<boost::shared_ptr<VirtualShape> > &getComponents();
  void updateComponents(double time);
  void fillGeometry();
  void fillEvent(const edm::Event& event);
  bool findTrajectory(const edm::Event& event, boost::shared_ptr<Track> track, int id);
  void findBoundaryT(minmax &m, double t);
  void findBoundaryP(minmax &m, double x, double y, double z);
  void makeSupportStructuresVisible(bool visible);
  void makeStrawsVisibleBeforeStart(bool visible);
  void useHitColors(bool hitcolors);
  void useTrackColors(bool trackcolors);
  minmax getHitsBoundary() {return _hitsMinmax;}
  minmax getTracksBoundary() {return _tracksMinmax;}
  int getNumberHits() {return _numberHits;}
};

}
#endif
