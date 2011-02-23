//
// Container class for all particle tracks. Tracks are displayed via the EventDisplayPolyLine3D class (inherited from ROOT's TPolyLine3D class). The displayed length of the track depends is time-dependent.
//
// $Id: Track.h,v 1.4 2011/02/23 00:29:27 ehrlich Exp $
// $Author: ehrlich $ 
// $Date: 2011/02/23 00:29:27 $
//
// Original author Ralf Ehrlich
//

#ifndef TRACK_H
#define TRACK_H

#include "dict_classes/EventDisplayPolyLine3D.h"
#include "VirtualShape.h"
#include <TPad.h>
#include <TMath.h>
#include <vector>

namespace mu2e_eventdisplay
{

class Track: public VirtualShape 
{
  Track();
  Track(const Track &);
  Track& operator=(const Track &);

  struct pt{double x,y,z,t;};
  std::vector<pt> _pVec;
  boost::shared_ptr<EventDisplayPolyLine3D> _line;
  bool _trajectory;
  int _particleId;

  public:
  Track(double x1, double y1, double z1, double t1, 
        double x2, double y2, double z2, double t2,
        int particleId,
        const TGeoManager *geomanager, TGeoVolume *topvolume,
        const TObject *mainframe, const boost::shared_ptr<ComponentInfo> info):
        VirtualShape(geomanager, topvolume, info, false)
  {
    pt p;
    p.x=x1;    
    p.y=y1;    
    p.z=z1;    
    p.t=t1;   
    _pVec.push_back(p); 
    p.x=x2;    
    p.y=y2;    
    p.z=z2;    
    p.t=t2;    
    _pVec.push_back(p); 
    setStartTime(t1); 
    setEndTime(t2); 
    _particleId=particleId;
    _trajectory=false;
    _line=boost::shared_ptr<EventDisplayPolyLine3D>(new EventDisplayPolyLine3D(mainframe, _info));
    _line->Draw();
    start();
  }

  void addTrajectoryPoint(double x, double y, double z, int i, int numberPoints)
  {
    if(!_trajectory) {_pVec.clear();}
    pt p;
    p.x=x; p.y=y; p.z=z;
    p.t=getStartTime()+i*(getEndTime()-getStartTime())/(numberPoints-1);
    _pVec.push_back(p);
    if(!_trajectory) {start(); _trajectory=true;}
  }

  virtual ~Track() 
  { //_line gets deleted automatically via boost
  }

  void start()
  {
    _line->SetPolyLine(0); //needed, e.g if animation is restarted
    _line->SetPoint(0,_pVec[0].x,_pVec[0].y,_pVec[0].z);
    _line->SetLineColor(getColor());
    _line->SetLineWidth(1);
  }

  int getParticleId() {return _particleId;}

  void toForeground()
  {
    gPad->RecursiveRemove(_line.get());
    _line->Draw();
  }

  void update(double time)
  {
    if(time<getStartTime()) return;
    for(unsigned int i=1; i<_pVec.size(); i++)
    {
      if(time>_pVec[i].t || isnan(time)) 
      {
        _line->SetPoint(i,_pVec[i].x,_pVec[i].y,_pVec[i].z);
        _line->SetLineColor(getColor());
      }
      if(time<=_pVec[i].t)
      {
        double scale=(time-_pVec[i-1].t)/(_pVec[i].t-_pVec[i-1].t);
        double x_tmp=_pVec[i-1].x+(_pVec[i].x-_pVec[i-1].x)*scale;
        double y_tmp=_pVec[i-1].y+(_pVec[i].y-_pVec[i-1].y)*scale;
        double z_tmp=_pVec[i-1].z+(_pVec[i].z-_pVec[i-1].z)*scale;
        _line->SetPoint(i,x_tmp,y_tmp,z_tmp);
        _line->SetLineColor(getColor());
        break;
      }
    }
  }

};

}
#endif
