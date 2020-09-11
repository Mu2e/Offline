//
// Class for all hexagon structures, e.g. crystals in disk calorimeter. The structure is displayed via EventDisplayGeoVolumePgon (inherited from TGeoVolume) which holds a TGeoPgon. In order to allow the user to right-click the structure and get a contect menu, there are additional lines drawn via the EventDisplayPolyLine3D class (inherited from ROOT's TPolyLine3D class).
//
//
// Original author Ralf Ehrlich
//

#ifndef EventDisplay_src_Hexagon_h
#define EventDisplay_src_Hexagon_h

#include "EventDisplay/src/VirtualShape.h"
#include "EventDisplay/src/dict_classes/EventDisplayGeoVolumePgon.h"
#include "EventDisplay/src/dict_classes/EventDisplayPolyLine3D.h"
#include <TPad.h>
#include <TMath.h>
#include <iostream>
#include <vector>

namespace mu2e_eventdisplay
{

class Hexagon: public VirtualShape
{
  Hexagon();
  Hexagon(const Hexagon &);
  Hexagon& operator=(const Hexagon &);

  //bare pointer needed, since ROOT manages this object
  EventDisplayGeoVolumePgon *_volume;
  TGeoRotation *_rotation;
  TGeoTranslation *_translation;
  std::vector<boost::shared_ptr<EventDisplayPolyLine3D> > _lines;

  struct points{double x,y,z;};

  public:
  Hexagon(double x, double y, double z, double rmax, double halflength, double phiOffset, double time,
           const TGeoManager *geomanager, TGeoVolume *topvolume,
           EventDisplayFrame *mainframe, const boost::shared_ptr<ComponentInfo> info,
           bool isGeometry) : 
           VirtualShape(geomanager, topvolume, mainframe, info, isGeometry)
  {
    setStartTime(time);
    _notDrawn=true;

    _volume = new EventDisplayGeoVolumePgon(phiOffset, rmax, halflength, mainframe, _info);
    _volume->SetVisibility(0);
    _volume->SetLineWidth(1);
    _translation = new TGeoTranslation(x,y,z);
    int i=0;
    if(_topvolume->GetNodes()) i=_topvolume->GetNodes()->GetEntries();
    _topvolume->AddNode(_volume, i, _translation);

    _lines.resize(8);
    for(int i=0; i<8; i++)
    {
      _lines[i] = boost::shared_ptr<EventDisplayPolyLine3D>(new EventDisplayPolyLine3D(mainframe, _info));
      _lines[i]->SetLineWidth(1);
    }
    for(int phi=phiOffset, i=0; phi<=phiOffset+360; phi+=60, i++)
    {
      double phiRad=phi*TMath::Pi()/180.0;
      double px=rmax*sin(phiRad);
      double py=rmax*cos(phiRad);
      _lines[0]->SetPoint(i,x,y,z-halflength);
      _lines[1]->SetPoint(i,x,y,z+halflength);
      if(i<6)
      {
        _lines[i+2]->SetPoint(0,x+px,y+py,z-halflength);
        _lines[i+2]->SetPoint(1,x+px,y+py,z+halflength);
      }
    }

    start();
  }

  virtual ~Hexagon()
  {
    _lines.clear();   //deletes all lines
    _volume->SetVisibility(0);  //can only be deleted by deleting TGeoManager
                                //deleting TGeoManager should also delete
                                //TGeoRotation, and TGeoTranslation
  }

  void start()
  {
    resetFilter();
    _volume->SetVisibility(0);
    if(_notDrawn==false)
    {
      for(int i=0; i<8; i++)
      {
        //gPad->RecursiveRemove(_lines[i].get());
        gPad->GetListOfPrimitives()->Remove(_lines[i].get());
      }
    }
    _notDrawn=true;
  }

  void toForeground()
  {
    if(_notDrawn==false)
    {
      for(int i=0; i<8; i++)
      {
        //gPad->RecursiveRemove(_lines[i].get());
        gPad->GetListOfPrimitives()->Remove(_lines[i].get());
        _lines[i]->Draw();
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

    _volume->SetVisibility(1);
    _volume->SetLineColor(getColor());
    _volume->SetFillColor(getColor());
    for(int i=0; i<8; i++)
    {
      _lines[i]->SetLineColor(getColor());
      if(_notDrawn) _lines[i]->Draw();
    }
    _notDrawn=false;
  }

  void makeGeometryVisible(bool visible)
  {
    if(visible)
    {
      _volume->SetVisibility(1);
      _volume->SetLineColor(getDefaultColor());
      _volume->SetFillColor(getDefaultColor());
      if(_notDrawn==true)
      {
        for(int i=0; i<8; i++)
        {
          _lines[i]->SetLineColor(getDefaultColor());
          _lines[i]->Draw();
        }
      }
      _notDrawn=false;
    }
    else start();
  }

  ClassDef(Hexagon,1);
};

}
#endif /* EventDisplay_src_Hexagon_h */
