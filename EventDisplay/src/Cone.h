//
// Class for all static (i.e. time-independent) cone structures. Now only used for ProtonAbsorber. The structure is displayed via EventDisplayGeoVolumeCone (inherited from TGeoVolume) which holds a TGeoCone.
//
//
// Original author MyeongJae Lee, based on Ralf Ehrlich's Tube.h.
//

#ifndef EventDisplay_src_Cone_h
#define EventDisplay_src_Cone_h

#include "Offline/EventDisplay/src/VirtualShape.h"
#include "Offline/EventDisplay/src/dict_classes/EventDisplayGeoVolumeCone.h"
#include "Offline/EventDisplay/src/dict_classes/EventDisplayPolyLine3D.h"
#include <TMath.h>
#include <iostream>
#include <map>

namespace mu2e_eventdisplay
{

class Cone: public VirtualShape
{
  Cone();
  Cone(const Cone &);
  Cone& operator=(const Cone &);

  //bare pointer needed, since ROOT manages this object
  EventDisplayGeoVolumeCone *_volume;
  TGeoRotation *_rotation;
  TGeoCombiTrans *_translation;
  struct line_struct
  {
    double x1, y1, z1, x2, y2, z2;
    boost::shared_ptr<EventDisplayPolyLine3D> line;
  };
  std::vector<line_struct> _lines;
  /*
  struct layersegment_struct
  {
    int layer, segment, direction;
    bool operator<(const layersegment_struct &rhs) const
    {
      if(this->layer<rhs.layer) return true;
      if(this->layer>rhs.layer) return false;
      if(this->segment<rhs.segment) return true;
      if(this->segment>rhs.segment) return false;
      if(this->direction<rhs.direction) return true;
      return false;
    }
    bool operator==(const layersegment_struct &rhs) const
    {
      if(this->layer==rhs.layer && this->segment==rhs.segment && this->direction==rhs.direction) return true;
      return false;
    }
  };
  std::map<layersegment_struct,line_struct> _lines;
*/
  void addline(double x1, double y1, double z1, double x2, double y2, double z2,
               int layer, int segment, int direction, EventDisplayFrame *mainframe)
  {
    line_struct newline;
    newline.x1=x1;
    newline.y1=y1;
    newline.z1=z1;
    newline.x2=x2;
    newline.y2=y2;
    newline.z2=z2;
    newline.line=boost::shared_ptr<EventDisplayPolyLine3D>(new EventDisplayPolyLine3D(mainframe, _info));
    newline.line->SetLineWidth(1);
    newline.line->SetPoint(0,x1,y1,z1);
    newline.line->SetPoint(1,x2,y2,z2);
/*    layersegment_struct layersegment;
    layersegment.layer=layer;
    layersegment.segment=segment;
    layersegment.direction=direction;
    _lines[layersegment]=newline;*/
    _lines.push_back(newline);
  }

  public:

  Cone(double x, double y, double z,
           double phi, double theta, double psi, double halflength,
           double innerRadius1, double outerRadius1,
           double innerRadius2, double outerRadius2,
           double startTime,
           const TGeoManager *geomanager, TGeoVolume *topvolume,
           EventDisplayFrame *mainframe, const boost::shared_ptr<ComponentInfo> info,
           bool isGeometry):
           VirtualShape(geomanager, topvolume, mainframe, info, isGeometry)
  {
    setStartTime(startTime);
    _notDrawn=true;

    _volume = new EventDisplayGeoVolumeCone(innerRadius1, outerRadius1, innerRadius2, outerRadius2, halflength, mainframe, _info);
    _volume->SetVisibility(0);
    _volume->SetLineWidth(1);
    _rotation = new TGeoRotation("",phi*180.0/TMath::Pi(),theta*180.0/TMath::Pi(),psi*180.0/TMath::Pi());
    _translation = new TGeoCombiTrans(x,y,z,_rotation);
    int i=0;
    if(_topvolume->GetNodes()) i=_topvolume->GetNodes()->GetEntries();
    _topvolume->AddNode(_volume, i, _translation);


    int nseg=_geomanager->GetNsegments();
    int nlayers=TMath::CeilNint((outerRadius1-innerRadius1)/500.0);
    double st=sin(theta);
    double ct=cos(theta);
    double sp=sin(phi);
    double cp=cos(phi);
    double ss=sin(psi);
    double cs=cos(psi);
    for(int i=0; i<=nlayers; i++)
    {
      double r1=i*(outerRadius1-innerRadius1)/nlayers+innerRadius1;
      double r2=i*(outerRadius2-innerRadius2)/nlayers+innerRadius2;
      for(int j=0; j<nseg; j++)
      {
        if(i==0 && innerRadius1==0) j=nseg;
        double azimuth=j*2*TMath::Pi()/nseg;

        //Start with an unrotated cylinder with its center at (0,0,0).
        //The vectors to points on the end planes of the cylinder
        //before rotation are (ex,ey,-halflength) and (ex,ey,halflength).
        double dx1 = cos(azimuth+TMath::Pi()/2.0)*r1;
        double dy1 = sin(azimuth+TMath::Pi()/2.0)*r1;
        double dx2 = cos(azimuth+TMath::Pi()/2.0)*r2;
        double dy2 = sin(azimuth+TMath::Pi()/2.0)*r2;
        double dz1=-halflength;
        double dz2=halflength;

        //After the rotation, the vectors from the center of the cylinder
        //to the points on the end planes of the cylinder
        //are (rx1,ry1,-halflength) and (rx2,ry2,halflength).
        double rx1,ry1,rz1;
        double rx2,ry2,rz2;
        rotate(dx1,dy1,dz1,  rx1,ry1,rz1,  sp,cp,st,ct,ss,cs);
        rotate(dx2,dy2,dz2,  rx2,ry2,rz2,  sp,cp,st,ct,ss,cs);

        //After the translation (i.e. when the center of the cylinder moves
        //from (0,0,0) to (x,y,z)), the points of the cylinder points move to
        //(x+rx1,y+ry1,z+rz1) and (x+rx2,y+ry2,z+rz2).
        addline(rx1+x, ry1+y, rz1+z,   rx2+x, ry2+y, rz2+z,   i, j, 0, mainframe);
      }
    }

/*
    for(int i=0; i<=nlayers; i++)
    {
      if(i==0 && innerRadius==0) continue;
      for(int j=0; j<nseg; j++)
      {
        layersegment_struct layersegment0;
        layersegment0.layer=i;
        layersegment0.segment=j;
        layersegment0.direction=0;
        const line_struct &referenceline0=_lines[layersegment0];

        layersegment_struct layersegment1;
        layersegment1.layer=i;
        layersegment1.segment=j+1;
        layersegment1.direction=0;
        if(layersegment1.segment==nseg) layersegment1.segment=0;
        const line_struct &referenceline1=_lines[layersegment1];
        addline(referenceline0.x1, referenceline0.y1, referenceline0.z1,
                referenceline1.x1, referenceline1.y1, referenceline1.z1,
                i, j, 1, mainframe);
        addline(referenceline0.x2, referenceline0.y2, referenceline0.z2,
                referenceline1.x2, referenceline1.y2, referenceline1.z2,
                i, j, 2, mainframe);

        if(i==0) continue;
        layersegment_struct layersegment2;
        layersegment2.layer=i-1;
        layersegment2.segment=j;
        layersegment2.direction=0;
        if(i==1 && innerRadius==0) layersegment2.segment=nseg;
        const line_struct &referenceline2=_lines[layersegment2];
        addline(referenceline0.x1, referenceline0.y1, referenceline0.z1,
                referenceline2.x1, referenceline2.y1, referenceline2.z1,
                i, j, 3, mainframe);
        addline(referenceline0.x2, referenceline0.y2, referenceline0.z2,
                referenceline2.x2, referenceline2.y2, referenceline2.z2,
                i, j, 4, mainframe);
      }
    }
*/

    start();
  }

  virtual ~Cone()
  {
    _lines.clear();   //deletes all lines
    _volume->SetVisibility(0);  //can only be deleted by deleting TGeoManager
                                //deleting TGeoManager should also delete
                                //TGeoRotation, and TGeoTranslation
  }

  void start()
  {
    //std::map<layersegment_struct,line_struct>::iterator iter;
    _volume->SetVisibility(0);
    if(_notDrawn==false)
    {
      std::vector<line_struct>::iterator iter;
      for(iter=_lines.begin(); iter!=_lines.end(); iter++)
      {
        //line_struct &l=iter->second;
        line_struct &l=*iter;
        //gPad->RecursiveRemove(l.line.get());
        gPad->GetListOfPrimitives()->Remove(l.line.get());
      }
    }
    _notDrawn=true;
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
    //std::map<layersegment_struct,line_struct>::iterator iter;
    std::vector<line_struct>::iterator iter;
    for(iter=_lines.begin(); iter!=_lines.end(); iter++)
    {
      //line_struct &l=iter->second;
      line_struct &l=*iter;
      l.line->SetLineColor(getColor());
      if(_notDrawn==true) l.line->Draw();
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
        std::vector<line_struct>::iterator iter;
        for(iter=_lines.begin(); iter!=_lines.end(); iter++)
        {
          //line_struct &l=iter->second;
          line_struct &l=*iter;
          l.line->SetLineColor(getDefaultColor());
          l.line->Draw();
        }
      }
      _notDrawn=false;
    }
    else start();
  }

  ClassDef(Cone,1);
};

}
#endif /* EventDisplay_src_Cone_h */
