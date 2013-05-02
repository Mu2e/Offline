//
// Container class for all detector straws. Straws are displayed via the EventDisplayPolyLine3D class (inherited from ROOT's TPolyLine3D class). Straws which are hit are drawn in a particular color which depends on the hit time.
//
// $Id: Straw.h,v 1.12 2013/05/02 06:03:41 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2013/05/02 06:03:41 $
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
  bool _notDrawn;
  int  _hitNumber;
  bool _invisible;

  public:

  Straw(double x, double y, double z, double t1,
        double theta, double phi, double halflength,
        const TGeoManager *geomanager, TGeoVolume *topvolume,
        EventDisplayFrame *mainframe, const boost::shared_ptr<ComponentInfo> info):
        VirtualShape(geomanager, topvolume, mainframe, info, true)
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
    _line->SetLineColor(getDefaultColor());
    _line->SetPoint(0,x1,y1,z1);
    _line->SetPoint(1,x2,y2,z2);
    start();
  }

  virtual ~Straw()
  {
  }

  void start()
  {
    _invisible=false;
    if(_notDrawn==false) gPad->RecursiveRemove(_line.get());
    _notDrawn=true;
  }

  void toForeground()
  {
    if(_notDrawn==false)
    {
      gPad->RecursiveRemove(_line.get());
      _line->Draw();
    }
  }

  void update(double time)
  {
    if(time<getStartTime() || isnan(getStartTime())) return;

    if(_invisible)
    {
      if(_notDrawn==false) gPad->RecursiveRemove(_line.get());
      _notDrawn=true;
      return;
    }

    _line->SetLineColor(getColor());
    if(_notDrawn)
    {
      if(_notDrawn==true) _line->Draw();
      _notDrawn=false;
    }
  }
 
  void setInvisible(bool invisible)
  {
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

};

}
#endif /* EventDisplay_src_Straw_h */
