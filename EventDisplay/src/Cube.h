//
// Class for all cube structures, e.g. vanes, crystals. The structure is displayed via EventDisplayGeoVolumeBox (inherited from TGeoVolume) which holds a TGeoBox. In order to allow the user to right-click the structure and get a contect menu, there are additional lines drawn via the EventDisplayPolyLine3D class (inherited from ROOT's TPolyLine3D class).
//
//
// Original author Ralf Ehrlich
//

#ifndef EventDisplay_src_Cube_h
#define EventDisplay_src_Cube_h

#include "Offline/EventDisplay/src/VirtualShape.h"
#include "Offline/EventDisplay/src/dict_classes/EventDisplayGeoVolumeBox.h"
#include "Offline/EventDisplay/src/dict_classes/EventDisplayPolyLine3D.h"
#include <TPad.h>
#include <TMath.h>
#include <iostream>
#include <vector>

namespace mu2e_eventdisplay
{

class Cube: public VirtualShape
{
  Cube();
  Cube(const Cube &);
  Cube& operator=(const Cube &);

  //bare pointer needed, since ROOT manages this object
  EventDisplayGeoVolumeBox *_volume;
  TGeoRotation *_rotation;
  TGeoCombiTrans *_translation;
  struct line_struct
  {
    boost::shared_ptr<EventDisplayPolyLine3D> line;
  };
  std::vector<line_struct> _lines;

  struct points{double x,y,z;};

  points rotateAndTranslate(double x, double y, double z,
                            double dx, double dy, double dz,
                            double phi, double theta, double psi)
  {
    double st=sin(theta);
    double ct=cos(theta);
    double sp=sin(phi);
    double cp=cos(phi);
    double ss=sin(psi);
    double cs=cos(psi);

    //Start with the cylinder centered at (0,0,0).
    //Before the rotation, the vector from the center of the cube
    //to its corner is (dx,dy,dz).
    //After the rotation, the vector from the center of the cube
    //to its corner is (rx,ry,rz).
    double rx,ry,rz;
    rotate(dx,dy,dz,  rx,ry,rz,  sp,cp,st,ct,ss,cs);

    //After the translation (i.e. when the center of the cube moves
    //from (0,0,0) to (x,y,z)), the points of the cube points move to
    //(x+rx,y+ry,z+rz).
    points to_return;
    to_return.x=x+rx;
    to_return.y=y+ry;
    to_return.z=z+rz;

    return to_return;
  }


  public:
  Cube(double x, double y, double z, double dx, double dy, double dz,
           double phi, double theta, double psi, double time,
           const TGeoManager *geomanager, TGeoVolume *topvolume,
           EventDisplayFrame *mainframe, const boost::shared_ptr<ComponentInfo> info,
           bool isGeometry) :
           VirtualShape(geomanager, topvolume, mainframe, info, isGeometry)
  {
    setStartTime(time);
    _notDrawn=true;

    _volume = new EventDisplayGeoVolumeBox(dx, dy, dz, mainframe, _info);
    _volume->SetVisibility(0);
    _volume->SetLineWidth(2);
    _rotation = new TGeoRotation("",phi*180.0/TMath::Pi(),theta*180.0/TMath::Pi(),psi*180.0/TMath::Pi());
    _translation = new TGeoCombiTrans(x,y,z,_rotation);
    int i=0;
    if(_topvolume->GetNodes()) i=_topvolume->GetNodes()->GetEntries();
    _topvolume->AddNode(_volume, i, _translation);

    //Start with the cylinder centered at (0,0,0).
    //Before the rotation, the vectors from the center of the cube
    //to its corners are
    //(dx,dy,dz), (-dx,dy,dz), (-dx,-dy,dz), (dx,-dy,dz), and so on.

    //Now rotate this cylinder with the angles phi, theta, psi,
    //and move the center of the cylinder from (0,0,0) to (x,y,z)
    points p[8];
    p[0]=rotateAndTranslate(x, y, z,  dx,  dy,  dz, phi, theta, psi);
    p[1]=rotateAndTranslate(x, y, z, -dx,  dy,  dz, phi, theta, psi);
    p[2]=rotateAndTranslate(x, y, z, -dx, -dy,  dz, phi, theta, psi);
    p[3]=rotateAndTranslate(x, y, z,  dx, -dy,  dz, phi, theta, psi);
    p[4]=rotateAndTranslate(x, y, z,  dx,  dy, -dz, phi, theta, psi);
    p[5]=rotateAndTranslate(x, y, z, -dx,  dy, -dz, phi, theta, psi);
    p[6]=rotateAndTranslate(x, y, z, -dx, -dy, -dz, phi, theta, psi);
    p[7]=rotateAndTranslate(x, y, z,  dx, -dy, -dz, phi, theta, psi);

    line_struct newline;
    for(int i=0; i<4; i++)
    {
      newline.line=boost::shared_ptr<EventDisplayPolyLine3D>(new EventDisplayPolyLine3D(mainframe, _info));
      newline.line->SetLineWidth(2);
      newline.line->SetPoint(0,p[i].x,p[i].y,p[i].z);
      newline.line->SetPoint(1,p[i+4].x,p[i+4].y,p[i+4].z);
      _lines.push_back(newline);
    }
    for(int i=0; i<4; i++)
    {
      newline.line=boost::shared_ptr<EventDisplayPolyLine3D>(new EventDisplayPolyLine3D(mainframe, _info));
      newline.line->SetLineWidth(2);
      newline.line->SetPoint(0,p[i].x,p[i].y,p[i].z);
      if(i<3) newline.line->SetPoint(1,p[i+1].x,p[i+1].y,p[i+1].z);
      else newline.line->SetPoint(1,p[0].x,p[0].y,p[0].z);
      _lines.push_back(newline);
    }
    for(int i=4; i<8; i++)
    {
      newline.line=boost::shared_ptr<EventDisplayPolyLine3D>(new EventDisplayPolyLine3D(mainframe, _info));
      newline.line->SetLineWidth(2);
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
    resetFilter();
    _volume->SetVisibility(0);
    if(_notDrawn==false)
    {
      std::vector<line_struct>::iterator iter;
      for(iter=_lines.begin(); iter!=_lines.end(); iter++)
      {
        line_struct &l=*iter;
//        gPad->RecursiveRemove(l.line.get());
        gPad->GetListOfPrimitives()->Remove(l.line.get());
      }
    }
    _notDrawn=true;
  }

  void toForeground()
  {
    if(_notDrawn==false)
    {
      std::vector<line_struct>::iterator iter;
      for(iter=_lines.begin(); iter!=_lines.end(); iter++)
      {
        line_struct &l=*iter;
//        gPad->RecursiveRemove(l.line.get());
        gPad->GetListOfPrimitives()->Remove(l.line.get());
        l.line->Draw();
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
    std::vector<line_struct>::iterator iter;
    for(iter=_lines.begin(); iter!=_lines.end(); iter++)
    {
      line_struct &l=*iter;
      l.line->SetLineColor(getColor());
      if(_notDrawn) l.line->Draw();
    }
    _notDrawn=false;
  }

  void makeGeometryVisible(bool visible)
  {
    if(visible)
    {
      std::vector<line_struct>::iterator iter;
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
    else start();
  }

  ClassDef(Cube,1);
};

}
#endif /* EventDisplay_src_Cube_h */
