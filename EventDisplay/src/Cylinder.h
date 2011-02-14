//
// Template class for all static (i.e. time-independent) cylinder structures, e.g. TTracker, Target. The structure is displayed via TGeoVolumeType (inherited from TGeoVolume) which holds a TGeoTube. In order to allow the user to right-click the structure and get a contect menu, there are additional lines drawn via the TPolyLine3DType class (inherited from ROOT's TPolyLine3D class). 
//
// $Id: Cylinder.h,v 1.3 2011/02/14 03:45:02 ehrlich Exp $
// $Author: ehrlich $ 
// $Date: 2011/02/14 03:45:02 $
//
// Original author Ralf Ehrlich
//

#ifndef CYLINDER_TEMPLATE_H
#define CYLINDER_TEMPLATE_H

#include "dict_classes/TGeoVolumeSupport.h"
#include "dict_classes/TGeoVolumeTarget.h"
#include "dict_classes/TPolyLine3DSupport.h"
#include "dict_classes/TPolyLine3DTarget.h"
#include <TMath.h>
#include "VirtualShape.h"
#include <iostream>
#include <vector>

namespace mu2e_eventdisplay
{

template<typename TGeoVolumeType, typename TPolyLine3DType>
class Cylinder: public VirtualShape 
{
  Cylinder();
  Cylinder(const Cylinder &);
  Cylinder& operator=(const Cylinder &);

  //bare pointer needed, since ROOT manages this object 
  TGeoVolumeType *_volume;
  TGeoRotation *_rotation;
  TGeoCombiTrans *_translation;
  struct line_struct 
  {
    double x1, y1, z1, x2, y2, z2;
    boost::shared_ptr<TPolyLine3DType> line;
  };
  std::vector<line_struct> _lines;
 
  public:

  Cylinder(double x, double y, double z,
           double phi, double theta, double psi, double halflength, 
           double innerRadius, double outerRadius, 
           const TGeoManager *geomanager, TGeoVolume *topvolume, 
           const TObject *mainframe, const boost::shared_ptr<ComponentInfo> info,
           bool defaultVisibility):
           VirtualShape(geomanager, topvolume, info, true)
  {
    setStartTime(NAN);
    setDefaultVisibility(defaultVisibility);

    _volume = new TGeoVolumeType(innerRadius, outerRadius, halflength, mainframe, _info);
    _volume->SetVisibility(0);
    _rotation = new TGeoRotation("",phi*180.0/TMath::Pi(),theta*180.0/TMath::Pi(),psi*180.0/TMath::Pi());
    _translation = new TGeoCombiTrans(x,y,z,_rotation);
    int i=0;
    if(_topvolume->GetNodes()) i=_topvolume->GetNodes()->GetEntries();
    _topvolume->AddNode(_volume, i, _translation);

    int nseg=_geomanager->GetNsegments();
    int nlayers=TMath::CeilNint((outerRadius-innerRadius)/200.0);
    double st=sin(theta);
    double ct=cos(theta);
    double sp=sin(phi);
    double cp=cos(phi);
    double ss=sin(psi);
    double cs=cos(psi);
    for(int i=0; i<=nlayers; i++)
    {
      double r=i*(outerRadius-innerRadius)/nlayers+innerRadius;
      for(int j=0; j<nseg; j++)
      {
        if(r==0) j=nseg;
        double azimuth=j*2*TMath::Pi()/nseg;

        //Start with an unrotated cylinder with its center at (0,0,0).
        //The vectors to points on the end planes of the cylinder 
        //before rotation are (ex,ey,-halflength) and (ex,ey,halflength).
        double dx = cos(azimuth+TMath::Pi()/2.0)*r;
        double dy = sin(azimuth+TMath::Pi()/2.0)*r;
        double dz1=-halflength;
        double dz2=halflength;

        //After the rotation, the vectors from the center of the cylinder
        //to the points on the end planes of the cylinder 
        //are (rx1,ry1,-halflength) and (rx2,ry2,halflength).
        double rx1,ry1,rz1;
        double rx2,ry2,rz2;
        rotate(dx,dy,dz1,  rx1,ry1,rz1,  sp,cp,st,ct,ss,cs);
        rotate(dx,dy,dz2,  rx2,ry2,rz2,  sp,cp,st,ct,ss,cs);

        //After the translation (i.e. when the center of the cylinder moves 
        //from (0,0,0) to (x,y,z)), the points of the cylinder points move to
        //(x+rx1,y+ry1,z+rz1) and (x+rx2,y+ry2,z+rz2).
        line_struct newline;
        newline.x1=rx1+x;
        newline.y1=ry1+y;
        newline.z1=rz1+z;
        newline.x2=rx2+x;
        newline.y2=ry2+y;
        newline.z2=rz2+z;
        newline.line=boost::shared_ptr<TPolyLine3DType>(new TPolyLine3DType(mainframe, _info));
        newline.line->Draw();
        _lines.push_back(newline);
      }
    }
    start();
  }

  virtual ~Cylinder() 
  {
    _lines.clear();   //deletes all lines
    _volume->SetVisibility(0);  //can only be deleted by deleting TGeoManager
                                //deleting TGeoManager should also delete 
                                //TGeoRotation, and TGeoTranslation
  }

  void start()
  {
    _volume->SetLineWidth(3);
    _volume->SetLineColor(getDefaultColor());
    _volume->SetFillColor(getDefaultColor());
    _volume->SetVisibility(0);
    if(getDefaultVisibility())
    {
      _volume->SetVisibility(1);
    }

    typedef std::vector<line_struct> vector_type;
    typedef typename vector_type::iterator iter_type;
    iter_type iter;
    for(iter=_lines.begin(); iter!=_lines.end(); iter++)
    {
      line_struct &l=*iter;
      l.line->SetPolyLine(0); //needed if animation is restarted
      l.line->SetPoint(0,l.x1,l.y1,l.z1);
      l.line->SetLineColor(getDefaultColor());
      l.line->SetLineWidth(3);
      if(getDefaultVisibility())
      {
        l.line->SetPoint(1,l.x2,l.y2,l.z2);
      }
    }
  }

  void update(double time)
  {//no animation here
  }
};

typedef Cylinder<TGeoVolumeSupport,TPolyLine3DSupport> SupportTTracker;
typedef Cylinder<TGeoVolumeTarget,TPolyLine3DTarget> Target;

}
#endif
