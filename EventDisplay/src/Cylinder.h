//
// Template class for all static (i.e. time-independent) cylinder structures, e.g. TTracker, Target. The structure is displayed via TGeoVolumeType (inherited from TGeoVolume) which holds a TGeoTube. In order to allow the user to right-click the structure and get a contect menu, there are additional lines drawn via the TPolyLine3DType class (inherited from ROOT's TPolyLine3D class). 
//
// $Id: Cylinder.h,v 1.1 2011/02/03 07:37:03 ehrlich Exp $
// $Author: ehrlich $ 
// $Date: 2011/02/03 07:37:03 $
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

  double _x, _y, _z;  
  double _theta, _phi, _halflength, _innerRadius, _outerRadius;
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
           double theta, double phi, double halflength, 
           double innerRadius, double outerRadius, 
           const TGeoManager *geomanager, TGeoVolume *topvolume, 
           const TObject *mainframe, const boost::shared_ptr<ComponentInfo> info,
           bool defaultVisibility):
           VirtualShape(geomanager, topvolume, info, true),
           _x(x),_y(y),_z(z),
           _theta(theta),_phi(phi),_halflength(halflength), 
           _innerRadius(innerRadius),_outerRadius(outerRadius)
  {
    setStartTime(NAN);
    setDefaultVisibility(defaultVisibility);

    _volume = new TGeoVolumeType(_innerRadius, _outerRadius, _halflength, mainframe, _info);
    _volume->SetVisibility(0);
    _rotation = new TGeoRotation("",_phi,_theta,0);
    _translation = new TGeoCombiTrans(_x,_y,_z,_rotation);
    int i=0;
    if(_topvolume->GetNodes()) i=_topvolume->GetNodes()->GetEntries();
    _topvolume->AddNode(_volume, i, _translation);

    int nseg=_geomanager->GetNsegments();
    int nlayers=TMath::CeilNint((_outerRadius-_innerRadius)/200.0);
    for(int i=0; i<=nlayers; i++)
    {
      double st=sin(_theta);
      double ct=cos(_theta);
      double sp=sin(_phi+TMath::Pi()/2.0);
      double cp=cos(_phi+TMath::Pi()/2.0);
      double r=i*(_outerRadius-_innerRadius)/nlayers+_innerRadius;
      for(int j=0; j<nseg; j++)
      {
        double rx=0;
        double ry=0;
        double rz=0;
        if(r>0)
        {
          double azimuth=j*2*TMath::Pi()/nseg;
          //before rotation:
          double ex = cos(azimuth+TMath::Pi()/2.0)*r;
          double ey = sin(azimuth+TMath::Pi()/2.0)*r;
          //after rotation:
          rx = cp*ex - ct*sp*ey;
          ry = sp*ex + ct*cp*ey;
          rz =         st*ey;
        }
        else
        {
          j=nseg;
        }
        //after translation:
        line_struct newline;
        newline.x1=rx+_x+_halflength*st*sp;
        newline.y1=ry+_y-_halflength*st*cp;
        newline.z1=rz+_z+_halflength*ct;
        newline.x2=rx+_x-_halflength*st*sp;
        newline.y2=ry+_y+_halflength*st*cp;
        newline.z2=rz+_z-_halflength*ct;
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
