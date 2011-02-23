//
// Container class for all detector straws. Straws are displayed via the EventDisplayPolyLine3D class (inherited from ROOT's TPolyLine3D class). Straws which are hit are drawn in a particular color which depends on the hit time.
//
// $Id: Straw.h,v 1.5 2011/02/23 00:29:27 ehrlich Exp $
// $Author: ehrlich $ 
// $Date: 2011/02/23 00:29:27 $
//
// Original author Ralf Ehrlich
//

#ifndef STRAW_H
#define STRAW_H

#include "dict_classes/EventDisplayPolyLine3D.h"
#include <TPad.h>
#include <TMath.h>
#include "VirtualShape.h"
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
 
  public:

  Straw(double x, double y, double z, double t1,
        double theta, double phi, double halflength, 
        const TGeoManager *geomanager, TGeoVolume *topvolume, 
        const TObject *mainframe, const boost::shared_ptr<ComponentInfo> info, 
        bool defaultVisibility):
        VirtualShape(geomanager, topvolume, info, true)
  {
    setStartTime(t1); 
    setDefaultVisibility(defaultVisibility);
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
    if(getDefaultVisibility())
    {
      _line->SetLineColor(getDefaultColor());
      if(_notDrawn==true) _line->Draw();
      _notDrawn=false;
    }
    else
    {
      if(_notDrawn==false) gPad->RecursiveRemove(_line.get());
      _notDrawn=true;
    }
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
    _line->SetLineColor(getColor());
    if(_notDrawn)
    {
      if(_notDrawn==true) _line->Draw();
      _notDrawn=false;
    }
  }

};

}
#endif
