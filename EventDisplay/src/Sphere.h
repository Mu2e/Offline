//
// Class for all non-static (i.e. time-dependent) sphere structures, e.g. drift radii. The structure is displayed via EventDisplayPolyLine3D (inherited from ROOT's TPolyLine3D) lines which allows the user to right-click the structure and get a contect menu.
//
// $Id: Sphere.h,v 1.4 2014/02/22 01:52:18 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2014/02/22 01:52:18 $
//
// Original author Ralf Ehrlich
//

#ifndef EventDisplay_src_Sphere_h
#define EventDisplay_src_Sphere_h

#include "EventDisplay/src/VirtualShape.h"
#include "EventDisplay/src/dict_classes/EventDisplayPolyLine3D.h"
#include <TMath.h>
#include <iostream>
#include <map>

namespace mu2e_eventdisplay
{

class Sphere: public VirtualShape
{
  Sphere();
  Sphere(const Sphere &);
  Sphere& operator=(const Sphere &);

  std::vector<boost::shared_ptr<EventDisplayPolyLine3D> > _pVec;
  double _x0, _y0, _z0, _r0;

  void drawSphere(double radius)
  {
    double phi=TMath::Pi()/2.0;
    double theta=0;
    double psi=0;
    double sinphi=sin(phi);
    double cosphi=cos(phi);
    double sinpsi=sin(psi);
    double cospsi=cos(psi);
    for(int i=0; i<10; i++, theta+=TMath::Pi()/10.0)
    {
      double sintheta=sin(theta);
      double costheta=cos(theta);
      boost::shared_ptr<EventDisplayPolyLine3D> p(new EventDisplayPolyLine3D(_mainframe, _info));
      _pVec.push_back(p);
      p->SetLineColor(getColor());
      p->SetLineWidth(1);
      double angle=0;
      double angleStep=2.0*TMath::Pi()/16;
      for(int i=0; i<=16; i++, angle+=angleStep)
      {
        double x=radius*sin(angle);
        double y=radius*cos(angle);
        double z=0;

        double x_rotated, y_rotated, z_rotated;
        rotate(x,y,z,  x_rotated,y_rotated,z_rotated,
               sinphi, cosphi,
               sintheta, costheta,
               sinpsi, cospsi);
        p->SetPoint(i,x_rotated+_x0,y_rotated+_y0,z_rotated+_z0);
      }
      p->Draw();
    }

    theta=-TMath::Pi()/2.0+TMath::Pi()/10.0;
    for(int b=-4; b<5; b++, theta+=TMath::Pi()/10.0)
    {
      double newRadius=radius*cos(theta);
      boost::shared_ptr<EventDisplayPolyLine3D> p(new EventDisplayPolyLine3D(_mainframe, _info));
      _pVec.push_back(p);
      p->SetLineColor(getColor());
      p->SetLineWidth(1);
      double angle=0;
      double angleStep=2.0*TMath::Pi()/16;
      for(int i=0; i<=16; i++, angle+=angleStep)
      {
        double x=newRadius*sin(angle);
        double y=radius*sin(theta);
        double z=newRadius*cos(angle);
        p->SetPoint(i,x+_x0,y+_y0,z+_z0);
      }
      p->Draw();
    }
  }

  public:

  Sphere(double x, double y, double z, double radius, double t1,
         const TGeoManager *geomanager, TGeoVolume *topvolume,
         EventDisplayFrame *mainframe, const boost::shared_ptr<ComponentInfo> info,
         bool isGeometry) : 
         VirtualShape(geomanager, topvolume, mainframe, info, isGeometry),
         _x0(x), _y0(y), _z0(z), _r0(radius)
  {
    _notDrawn=true;
    setStartTime(t1);
  }

  ~Sphere()
  {
    removeSphere();
  }

  void start()
  {
    removeSphere();
  }

  void removeSphere()
  {
    _pVec.clear(); //should call the destructor of all EventDisplayPolyLine3Ds which are stored in the vector
                   //which calls the destructors of all TPolyLine3Ds
    _notDrawn=true;
  }

  void toForeground()
  {
    std::vector<boost::shared_ptr<EventDisplayPolyLine3D> >::iterator iter;
    for(iter=_pVec.begin(); iter!=_pVec.end(); iter++)
    {
      if(*iter)
      {
        //gPad->RecursiveRemove((*iter).get());
        gPad->GetListOfPrimitives()->Remove((*iter).get());
        (*iter)->Draw();
      }
    }
  }

  void update(double time)
  {
    if(time<getStartTime() || std::isnan(getStartTime()) || getStartTime()<_minTime || getStartTime()>_maxTime)
    {
      start();
      return;
    }

    if(time>=getStartTime())
    {
      if(_notDrawn)
      {
        drawSphere(_r0);
        _notDrawn=false;
      }
    }
  }

  ClassDef(Sphere,1);
};

}
#endif /* EventDisplay_src_Sphere_h */
