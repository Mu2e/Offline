//
// Virtual base class for all shapes.
// Container class for the geometry object(s) with information on how they are to be displayed and updated for specific times.
//
//
// Original author Ralf Ehrlich
//

#ifndef EventDisplay_src_VirtualShape_h
#define EventDisplay_src_VirtualShape_h

#include "boost/shared_ptr.hpp"
#include "EventDisplay/src/EventDisplayFrame.h"
#include "EventDisplay/src/dict_classes/ComponentInfo.h"
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoVolume.h>
#include <iostream>
#include <math.h>
#include <string>

namespace mu2e_eventdisplay
{

class VirtualShape : public TObject
{
  VirtualShape();
  VirtualShape(const VirtualShape &);
  VirtualShape& operator=(const VirtualShape &);

  double _startTime;
  double _endTime;
  bool   _isGeometry;
  int    _color;

  protected:
  const TGeoManager *_geomanager;
  TGeoVolume        *_topvolume;
  EventDisplayFrame *_mainframe;
  boost::shared_ptr<ComponentInfo> _info;
  bool   _notDrawn;
  double _minTime, _maxTime; //filter

  public:

  VirtualShape(const TGeoManager *geomanager, TGeoVolume *topvolume,
               EventDisplayFrame *mainframe,
               const boost::shared_ptr<ComponentInfo> info, bool isGeometry):
               _startTime(NAN),_endTime(NAN),
               _isGeometry(isGeometry), _color(kGray),
               _geomanager(geomanager), _topvolume(topvolume), 
               _mainframe(mainframe), _info(info)
  {}

  virtual ~VirtualShape()
  {
//TODO
  }

  boost::shared_ptr<ComponentInfo> getComponentInfo() const {return _info;}

  void setStartTime(double t) {_startTime=t;}
  void setEndTime(double t) {_endTime=t;}
  double getStartTime() const {return _startTime;}
  double getEndTime() const {return _endTime;}

  void setIsPermanent(bool g) {_isGeometry=g;}
  bool isGeometry() const {return _isGeometry;}

  void setColor(int c) {_color=c;}
  int  getColor() const {return _color;}
  int  getDefaultColor() const {return kGray;}

  virtual void start()=0;
  virtual void update(double time)=0;
  virtual void makeGeometryVisible(bool visible) {;}
  virtual void toForeground() {;}

  void resetFilter()
  {
    _minTime=NAN;
    _maxTime=NAN;
  }

  void setFilter(double minTime, double maxTime)
  {
    _minTime=minTime;
    _maxTime=maxTime;
  }

  static void rotate(double oldX, double oldY, double oldZ,
                     double &newX, double &newY, double &newZ,
                     double sinphi, double cosphi,
                     double sintheta, double costheta,
                     double sinpsi, double cospsi)
  {
    double &sp=sinphi;
    double &cp=cosphi;
    double &st=sintheta;
    double &ct=costheta;
    double &ss=sinpsi;
    double &cs=cospsi;
    newX = cs*cp*oldX-ct*sp*ss*oldX   -  ss*cp*oldY-ct*sp*cs*oldY  +  st*sp*oldZ;
    newY = cs*sp*oldX+ct*cp*ss*oldX   -  ss*sp*oldY+ct*cp*cs*oldY  -  st*cp*oldZ;
    newZ = ss*st*oldX                 +  cs*st*oldY                +     ct*oldZ;
  }
};

}
#endif /* EventDisplay_src_VirtualShape_h */
