//
// Container class for all detector straws. Straws are displayed via the TPolyLine3DStraw class (inherited from ROOT's TPolyLine3D class). Straws which are hit are drawn in a particular color which depends on the hit time.
//
// $Id: Straw.h,v 1.2 2011/01/29 02:14:20 ehrlich Exp $
// $Author: ehrlich $ 
// $Date: 2011/01/29 02:14:20 $
//
// Original author Ralf Ehrlich
//

#ifndef STRAW_H
#define STRAW_H

#include "dict_classes/TPolyLine3DStraw.h"
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

  double _x, _y, _z;  
  double _x2, _y2, _z2;  
  double _theta, _phi, _halflength;
  boost::shared_ptr<TPolyLine3DStraw> _line;
  bool _notDrawn;
 
  public:

  Straw(double x, double y, double z, double t1,
        double theta, double phi, double halflength, 
        const TGeoManager *geomanager, TGeoVolume *topvolume, 
        const TObject *mainframe, const boost::shared_ptr<ComponentInfo> info, 
        bool defaultVisibility):
        VirtualShape(geomanager, topvolume, info, true),
        _x(x),_y(y),_z(z),
        _theta(theta),_phi(phi),_halflength(halflength)
  {
    setStartTime(t1); 
    setDefaultVisibility(defaultVisibility);
    _notDrawn=true;
    _line=boost::shared_ptr<TPolyLine3DStraw>(new TPolyLine3DStraw(mainframe, _info));
    start();
  }

  virtual ~Straw() 
  {
  }

  void start()
  {
    double st=sin(_theta);
    double ct=cos(_theta);
    double sp=sin(_phi+TMath::Pi()/2.0);
    double cp=cos(_phi+TMath::Pi()/2.0);
    double x1=_x+_halflength*st*sp;
    double y1=_y-_halflength*st*cp;
    double z1=_z+_halflength*ct;
    _x2=_x-_halflength*st*sp;
    _y2=_y+_halflength*st*cp;
    _z2=_z-_halflength*ct;
    _line->SetPolyLine(0); //needed if animation is restarted
    _line->SetPoint(0,x1,y1,z1);
    _line->SetLineWidth(1);
    if(getDefaultVisibility())
    {
      _line->SetLineColor(getDefaultColor());
      _line->SetPoint(1,_x2,_y2,_z2);
    }

    if(!getDefaultVisibility() && isnan(getStartTime()))
    {
      if(_notDrawn==false) gPad->RecursiveRemove(_line.get());
      _notDrawn=true;
    }
    else
    {
      if(_notDrawn==true) _line->Draw();
      _notDrawn=false;
    }
  }

  void update(double time)
  {
    if(_notDrawn) return;
    if(time<getStartTime() || isnan(getStartTime())) return;

    _line->SetLineColor(getColor());
    _line->SetPoint(1,_x2,_y2,_z2);  //TODO: don't want to do this everytime
  }

};

}
#endif
