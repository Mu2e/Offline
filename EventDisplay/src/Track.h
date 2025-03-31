//
// Container class for all particle tracks. Tracks are displayed via the EventDisplayPolyLine3D class (inherited from ROOT's TPolyLine3D class). The displayed length of the track depends is time-dependent.
//
//
// Original author Ralf Ehrlich
//

#ifndef EventDisplay_src_Track_h
#define EventDisplay_src_Track_h

#include "Offline/EventDisplay/src/VirtualShape.h"
#include "Offline/EventDisplay/src/dict_classes/EventDisplayPolyLine3D.h"
#include <TMath.h>
#include <TPad.h>
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
  int _trackClass, _trackClassIndex;
  double _momentum;

  unsigned int _minPoints;
  double _minMomentum;
  bool _showElectrons;
  bool _showMuons;
  bool _showGammas;
  bool _showNeutrinos;
  bool _showNeutrons;
  bool _showOthers;

  public:
//constructor for tracks with start and end point only
//(detailed trajectory can be added later)
  Track(double x1, double y1, double z1, double t1,
        double x2, double y2, double z2, double t2,
        int particleId, int trackClass, int trackClassIndex,
        double momentum,
        const TGeoManager *geomanager, TGeoVolume *topvolume,
        EventDisplayFrame *mainframe, const boost::shared_ptr<ComponentInfo> info,
        bool isGeometry):
        VirtualShape(geomanager, topvolume, mainframe, info, isGeometry)
  {
    resetFilter();
    _momentum=momentum;
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
    _trackClass=trackClass;
    _trackClassIndex=trackClassIndex;
    _trajectory=false;
    _line=boost::shared_ptr<EventDisplayPolyLine3D>(new EventDisplayPolyLine3D(mainframe, _info));
    _line->SetLineWidth(3);
    _line->Draw();
    start();
  }
/*
//constructor for a track without any points
//trajectory points can be added later
  Track(int particleId, int trackClass, int trackClassIndex,
        double momentum,
        const TGeoManager *geomanager, TGeoVolume *topvolume,
        EventDisplayFrame *mainframe, const boost::shared_ptr<ComponentInfo> info):
        VirtualShape(geomanager, topvolume, mainframe, info, false)
  {
    resetFilter();
    _momentum=momentum;
    _particleId=particleId;
    _trackClass=trackClass;
    _trackClassIndex=trackClassIndex;
    _trajectory=false;
    _line=boost::shared_ptr<EventDisplayPolyLine3D>(new EventDisplayPolyLine3D(mainframe, _info));
    _line->Draw();
  }
*/
  void addTrajectoryPoint(double x, double y, double z, double t)
  {
    if(!_trajectory) {_pVec.clear();}
    pt p;
    p.x=x; p.y=y; p.z=z; p.t=t;
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
  }

  int getParticleId() {return _particleId;}
  int getTrackClass() {return _trackClass;}
  int getTrackClassIndex() {return _trackClassIndex;}

  void toForeground()
  {
    //gPad->RecursiveRemove(_line.get());
    gPad->GetListOfPrimitives()->Remove(_line.get());
    _line->Draw();
  }

  void resetFilter()
  {
    _minPoints=0;
    _minTime=NAN;
    _maxTime=NAN;
    _minMomentum=0;
    _showElectrons=true;
    _showMuons=true;
    _showGammas=true;
    _showNeutrinos=true;
    _showNeutrons=true;
    _showOthers=true;
  }

  void setFilter(unsigned int minPoints, double minTime, double maxTime, double minMomentum,
                 bool showElectrons, bool showMuons, bool showGammas, bool showNeutrinos, bool showNeutrons, bool showOthers)
  {
    _minPoints=minPoints;
    _minTime=minTime;
    _maxTime=maxTime;
    _minMomentum=minMomentum;
    _showElectrons=showElectrons;
    _showMuons=showMuons;
    _showGammas=showGammas;
    _showNeutrinos=showNeutrinos;
    _showNeutrons=showNeutrons;
    _showOthers=showOthers;
  }

  void update(double time)
  {
    //Filter
    bool toDelete=false;
    if(_pVec.size()<_minPoints) toDelete=true;
    if(getEndTime()<_minTime) toDelete=true;
    if(getStartTime()>_maxTime) toDelete=true;
    if(_momentum<_minMomentum) toDelete=true;
    switch(_particleId)
    {
                case   11:
                case  -11: if(!_showElectrons) toDelete=true; break;   //e+,e-
                case   13:
                case  -13: if(!_showMuons) toDelete=true; break;   //mu+,mu-
                case   22: if(!_showGammas) toDelete=true; break;   //gamma
                case 2112: if(!_showNeutrons) toDelete=true; break;   //n0
                case   12:
                case  -12:
                case   14:
                case  -14:
                case   16:
                case  -16: if(!_showNeutrinos) toDelete=true; break;   //neutrinos
                default  : if(!_showOthers) toDelete=true;
    };
    if(toDelete)
    {
      start();  //shortens the track to a length of 0
      return;
    }

    if(time<getStartTime())
    {
      start();  //shortens the track to a length of 0
      return;
    }

    for(unsigned int i=1; i<_pVec.size(); i++)
    {
      _line->SetLineColor(getColor());
      if(time>_pVec[i].t || std::isnan(time))
      {
        _line->SetPoint(i,_pVec[i].x,_pVec[i].y,_pVec[i].z);
      }
      if(time<=_pVec[i].t)
      {
        double scale=(time-_pVec[i-1].t)/(_pVec[i].t-_pVec[i-1].t);
        double x_tmp=_pVec[i-1].x+(_pVec[i].x-_pVec[i-1].x)*scale;
        double y_tmp=_pVec[i-1].y+(_pVec[i].y-_pVec[i-1].y)*scale;
        double z_tmp=_pVec[i-1].z+(_pVec[i].z-_pVec[i-1].z)*scale;
        _line->SetPoint(i,x_tmp,y_tmp,z_tmp);

        for(int j=i+1; j<_line->Size(); j++)
        {                 //needed, e.g if animation is restarted at an earlier time, which leaves
          _line->SetPoint(j,x_tmp,y_tmp,z_tmp);   //unused points which need to be overwritten
        }                                         //with the same point as the last one
        break;
      }
    }
  }

  ClassDef(Track,1);
};

}
#endif /* EventDisplay_src_Track_h */
