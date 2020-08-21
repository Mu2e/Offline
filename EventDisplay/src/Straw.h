//
// Container class for all detector straws. Straws are displayed via the EventDisplayPolyLine3D class (inherited from ROOT's TPolyLine3D class). Straws which are hit are drawn in a particular color which depends on the hit time.
//
//
// Original author Ralf Ehrlich
//

#ifndef EventDisplay_src_Straw_h
#define EventDisplay_src_Straw_h

#include "EventDisplay/src/VirtualShape.h"
#include "EventDisplay/src/dict_classes/EventDisplayPolyLine3D.h"
#include <TMath.h>
#include <TPad.h>
#include <iostream>
#include <vector>

namespace mu2e_eventdisplay
{

class Straw: public VirtualShape
{
  Straw();
  Straw(const Straw &);
  Straw& operator=(const Straw &);

  boost::shared_ptr<EventDisplayPolyLine3D> _line;
  int   _hitNumber;
  bool  _invisible;

  public:

  Straw(double x, double y, double z, double t1,
        double theta, double phi, double halflength,
        const TGeoManager *geomanager, TGeoVolume *topvolume,
        EventDisplayFrame *mainframe, const boost::shared_ptr<ComponentInfo> info,
        bool isGeometry):
        VirtualShape(geomanager, topvolume, mainframe, info, isGeometry)
  {
    setStartTime(t1);
    _hitNumber=-1;
    _notDrawn=true;
    _line=boost::shared_ptr<EventDisplayPolyLine3D>(new EventDisplayPolyLine3D(mainframe, _info));
    double st=sin(theta);
    double ct=cos(theta);
    double sp=sin(phi+TMath::Pi()/2.0);
    double cp=cos(phi+TMath::Pi()/2.0);
    double x1=x+halflength*st*sp;
    double y1=y-halflength*st*cp;
    double z1=z+halflength*ct;
    double x2=x-halflength*st*sp;
    double y2=y+halflength*st*cp;
    double z2=z-halflength*ct;
    _line->SetLineWidth(1);
    _line->SetPoint(0,x1,y1,z1);
    _line->SetPoint(1,x2,y2,z2);
    start();
  }

  virtual ~Straw()
  {    //since _line is a shared pointer, it gets automatically removed if Straw gets removed
  }

  void start()
  {
    resetFilter();
//    if(_notDrawn==false) //gPad->RecursiveRemove(_line.get());
    if(_notDrawn==false) gPad->GetListOfPrimitives()->Remove(_line.get());
    _notDrawn=true;
  }

  void toForeground()
  {
    if(_notDrawn==false)
    {
      //gPad->RecursiveRemove(_line.get());
      gPad->GetListOfPrimitives()->Remove(_line.get());
      _line->Draw();
    }
  }

  void update(double time)
  {
    if(time<getStartTime() || std::isnan(getStartTime()) || getStartTime()<_minTime || getStartTime()>_maxTime)
    {
      start();
      return;
    }

    _line->SetLineColor(getColor());
    if(_notDrawn)
    {
      if(_notDrawn==true) _line->Draw();
      _notDrawn=false;
    }
  }
 
  void resetFilter()
  {
    _minTime=NAN;
    _maxTime=NAN;
    _invisible=false;
  }

  void setFilter(double minTime, double maxTime, bool invisible)
  {
    _minTime=minTime;
    _maxTime=maxTime;
    _invisible=invisible; //true, if the flag value is not satisfied for this hit
  }
 
  void setHitNumber(int hitNumber)
  {
    _hitNumber=hitNumber;
  }

  int getHitNumber()
  {
    return _hitNumber;
  }

  ClassDef(Straw,1);
};

}
#endif /* EventDisplay_src_Straw_h */
