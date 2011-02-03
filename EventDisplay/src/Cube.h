//
// Template class for all cube structures, e.g. Vane, Crystal. The structure is displayed via TGeoVolumeType (inherited from TGeoVolume) which holds a TGeoBox. In order to allow the user to right-click the structure and get a contect menu, there are additional lines drawn via the TPolyLine3DType class (inherited from ROOT's TPolyLine3D class). 
//
// $Id: Cube.h,v 1.1 2011/02/03 07:37:03 ehrlich Exp $
// $Author: ehrlich $ 
// $Date: 2011/02/03 07:37:03 $
//
// Original author Ralf Ehrlich
//

#ifndef CUBE_TEMPLATE_H
#define CUBE_TEMPLATE_H

#include "dict_classes/TGeoVolumeCrystal.h"
#include "dict_classes/TGeoVolumeVane.h"
#include "dict_classes/TPolyLine3DCrystal.h"
#include "dict_classes/TPolyLine3DVane.h"
#include <TPad.h>
#include <TMath.h>
#include "VirtualShape.h"
#include <iostream>
#include <vector>

namespace mu2e_eventdisplay
{

template<typename TGeoVolumeType, typename TPolyLine3DType>
class Cube: public VirtualShape 
{
  Cube();
  Cube(const Cube &);
  Cube& operator=(const Cube &);

  //bare pointer needed, since ROOT manages this object 
  TGeoVolumeType *_volume;
  TGeoRotation *_rotation;
  TGeoCombiTrans *_translation;
  struct line_struct 
  {
    boost::shared_ptr<TPolyLine3DType> line;
  };
  std::vector<line_struct> _lines;
  bool _notDrawn;
 
  public:
  Cube(double x, double y, double z, double dx, double dy, double dz,
           double phi, double theta, double psi,
           double time, int color, 
           const TGeoManager *geomanager, TGeoVolume *topvolume, 
           const TObject *mainframe, const boost::shared_ptr<ComponentInfo> info,
           bool defaultVisibility):
           VirtualShape(geomanager, topvolume, info, true)
  {
    setStartTime(time);
    setColor(color);
    setDefaultVisibility(defaultVisibility);
    _notDrawn=true;

    _volume = new TGeoVolumeType(dx, dy, dz, mainframe, _info);
    _volume->SetVisibility(0);
    _volume->SetLineWidth(1);
    _volume->SetLineColor(getDefaultColor());
    _volume->SetFillColor(getDefaultColor());
    _rotation = new TGeoRotation("",phi*180.0/TMath::Pi(),theta*180.0/TMath::Pi(),psi*180.0/TMath::Pi());
    _translation = new TGeoCombiTrans(x,y,z,_rotation);
    int i=0;
    if(_topvolume->GetNodes()) i=_topvolume->GetNodes()->GetEntries();
    _topvolume->AddNode(_volume, i, _translation);

    double st=sin(theta);
    double ct=cos(theta);
    double sp=sin(phi);
    double cp=cos(phi);
    double ss=sin(psi);
    double cs=cos(psi);

    //vector from center of vane to corner after roatation:
    double rx = cs*cp*dx-ct*sp*ss*dx   -  ss*cp*dy-ct*sp*cs*dy  +  st*sp*dz;
    double ry = cs*sp*dx+ct*cp*ss*dx   -  ss*sp*dy+ct*cp*cs*dy  -  st*cp*dz;
    double rz = ss*st*dx               +  cs*st*dy              +     ct*dz;

    //after translation:
    struct points{double x,y,z;};
    points p[8];
    p[0].x=x+rx;
    p[0].y=y+ry;
    p[0].z=z+rz;
    p[1].x=x-rx;
    p[1].y=y+ry;
    p[1].z=z+rz;
    p[2].x=x-rx;
    p[2].y=y-ry;
    p[2].z=z+rz;
    p[3].x=x+rx;
    p[3].y=y-ry;
    p[3].z=z+rz;
    p[4].x=x+rx;
    p[4].y=y+ry;
    p[4].z=z-rz;
    p[5].x=x-rx;
    p[5].y=y+ry;
    p[5].z=z-rz;
    p[6].x=x-rx;
    p[6].y=y-ry;
    p[6].z=z-rz;
    p[7].x=x+rx;
    p[7].y=y-ry;
    p[7].z=z-rz;

    line_struct newline;
    for(int i=0; i<4; i++)
    {
      newline.line=boost::shared_ptr<TPolyLine3DType>(new TPolyLine3DType(mainframe, _info));
      newline.line->SetLineWidth(1);
      newline.line->SetLineColor(getDefaultColor());
      newline.line->SetPoint(0,p[i].x,p[i].y,p[i].z);
      newline.line->SetPoint(1,p[i+4].x,p[i+4].y,p[i+4].z);
      _lines.push_back(newline);
    }
    for(int i=0; i<4; i++)
    {
      newline.line=boost::shared_ptr<TPolyLine3DType>(new TPolyLine3DType(mainframe, _info));
      newline.line->SetLineWidth(1);
      newline.line->SetLineColor(getDefaultColor());
      newline.line->SetPoint(0,p[i].x,p[i].y,p[i].z);
      if(i<3) newline.line->SetPoint(1,p[i+1].x,p[i+1].y,p[i+1].z);
      else newline.line->SetPoint(1,p[0].x,p[0].y,p[0].z);
      _lines.push_back(newline);
    }
    for(int i=4; i<8; i++)
    {
      newline.line=boost::shared_ptr<TPolyLine3DType>(new TPolyLine3DType(mainframe, _info));
      newline.line->SetLineWidth(1);
      newline.line->SetLineColor(getDefaultColor());
      newline.line->SetPoint(0,p[i].x,p[i].y,p[i].z);
      if(i<7) newline.line->SetPoint(1,p[i+1].x,p[i+1].y,p[i+1].z);
      else newline.line->SetPoint(1,p[4].x,p[4].y,p[4].z);
      _lines.push_back(newline);
    }

    start();
  }

  virtual ~Cube() 
  {
    _lines.clear();   //deletes all lines
    _volume->SetVisibility(0);  //can only be deleted by deleting TGeoManager
                                //deleting TGeoManager should also delete 
                                //TGeoRotation, and TGeoTranslation
  }

  void start()
  {
    typedef std::vector<line_struct> vector_type;
    typedef typename vector_type::iterator iter_type;
    iter_type iter;
    if(getDefaultVisibility())
    {
      _volume->SetVisibility(1);
      _volume->SetLineColor(getDefaultColor());
      _volume->SetFillColor(getDefaultColor());
      if(_notDrawn==true)
      {
        for(iter=_lines.begin(); iter!=_lines.end(); iter++)
        {
          line_struct &l=*iter;
          l.line->SetLineColor(getDefaultColor());
          l.line->Draw();
        }
      }
      _notDrawn=false;
    }
    else
    {
      _volume->SetVisibility(0);
      if(_notDrawn==false)
      {
        for(iter=_lines.begin(); iter!=_lines.end(); iter++)
        {
          line_struct &l=*iter;
          gPad->RecursiveRemove(l.line.get());
        }
      }
      _notDrawn=true;
    }
  }

  void update(double time)
  {
    if(time<getStartTime() || isnan(getStartTime())) return;
    typedef std::vector<line_struct> vector_type;
    typedef typename vector_type::iterator iter_type;
    iter_type iter;
    for(iter=_lines.begin(); iter!=_lines.end(); iter++)
    {
      line_struct &l=*iter;
      l.line->SetLineColor(getColor());
      _volume->SetLineColor(getColor());
      _volume->SetFillColor(getColor());
    }
    if(_notDrawn)
    {
      _volume->SetVisibility(1);
      for(iter=_lines.begin(); iter!=_lines.end(); iter++)
      {
        line_struct &l=*iter;
        l.line->Draw();
      }
      _notDrawn=false;
    }
  }
};

typedef Cube<TGeoVolumeCrystal,TPolyLine3DCrystal> Crystal;
typedef Cube<TGeoVolumeVane,TPolyLine3DVane> Vane;

}
#endif
