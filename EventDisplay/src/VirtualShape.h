//
// Virtual base class for all shapes.
// Container class for the geometry object(s) with information on how they are to be displayed and updated for specific times.
//
// $Id: VirtualShape.h,v 1.4 2011/02/14 03:45:03 ehrlich Exp $
// $Author: ehrlich $ 
// $Date: 2011/02/14 03:45:03 $
//
// Original author Ralf Ehrlich
//

#ifndef VIRTUALSHAPE_H
#define VIRTUALSHAPE_H

#include <string>
#include <TGeoMatrix.h>
#include <TGeoVolume.h>
#include <TGeoManager.h>
#include <math.h>
#include <iostream>
#include "dict_classes/ComponentInfo.h"
#include "boost/shared_ptr.hpp"

namespace mu2e_eventdisplay
{

class VirtualShape : public TObject 
{
  VirtualShape();
  VirtualShape(const VirtualShape &);
  VirtualShape& operator=(const VirtualShape &);

  double _startTime;
  double _endTime;
  bool   _defaultVisibility;   //i.e. visibility of straws before hit, unhit straws, 
                              //permanent structures, etc.
  bool   _isGeometry;
  int    _color;

  protected:
  const TGeoManager *_geomanager;
  TGeoVolume        *_topvolume;
  boost::shared_ptr<ComponentInfo> _info;

  public:

  VirtualShape(const TGeoManager *geomanager, TGeoVolume *topvolume, 
               const boost::shared_ptr<ComponentInfo> info, bool isGeometry):
               _startTime(NAN),_endTime(NAN),
               _defaultVisibility(true), _isGeometry(isGeometry), _color(kGray),
               _geomanager(geomanager), _topvolume(topvolume), _info(info)
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

  void setDefaultVisibility(bool v) {_defaultVisibility=v;}
  void setIsPermanent(bool g) {_isGeometry=g;}
  bool getDefaultVisibility() const {return _defaultVisibility;}
  bool isGeometry() const {return _isGeometry;}

  void setColor(int c) {_color=c;}
  int  getColor() const {return _color;}
  int  getDefaultColor() const {return kGray;}

  virtual void start()=0;
  virtual void update(double time)=0;

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
    double &cs=sinpsi;
    newX = cs*cp*oldX-ct*sp*ss*oldX   -  ss*cp*oldY-ct*sp*cs*oldY  +  st*sp*oldZ;
    newY = cs*sp*oldX+ct*cp*ss*oldX   -  ss*sp*oldY+ct*cp*cs*oldY  -  st*cp*oldZ;
    newZ = ss*st*oldX                 +  cs*st*oldY                +     ct*oldZ;
  }
};

}
#endif
