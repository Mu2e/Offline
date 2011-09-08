//
// Class which sets up the 3D view of the main pad. It also provides functions that can handle user commands for zuum, rotate, etc.
//
// $Id: EventDisplayViewSetup.h,v 1.2 2011/09/08 03:54:45 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2011/09/08 03:54:45 $
//
// Original author Ralf Ehrlich
//

#ifndef EventDisplay_src_dict_classes_EventDisplayViewSetup_h
#define EventDisplay_src_dict_classes_EventDisplayViewSetup_h

#include <iostream>
#include <TPad.h>
#include <TView3D.h>
#include <TAxis3D.h>

namespace mu2e_eventdisplay
{

class EventDisplayViewSetup
{
  public:

  static void setup()
  {
    int irep=0;
    gPad->GetView()->SetView(180,60,90,irep);
    gPad->SetPhi(-90-180);
    gPad->SetTheta(90-60);
    gPad->GetView()->ShowAxis();
    TAxis3D::GetPadAxis(gPad)->SetLabelSize(0.025);
    TAxis3D::GetPadAxis(gPad)->SetTitleOffset(-0.5);
    TAxis3D::GetPadAxis(gPad)->SetXTitle("x [mm]");
    TAxis3D::GetPadAxis(gPad)->SetYTitle("y [mm]");
    TAxis3D::GetPadAxis(gPad)->SetZTitle("z [mm]");
    TAxis3D::GetPadAxis(gPad)->GetXaxis()->SetTitleColor(kRed);
    TAxis3D::GetPadAxis(gPad)->GetYaxis()->SetTitleColor(kGreen);
    TAxis3D::GetPadAxis(gPad)->GetZaxis()->SetTitleColor(kBlue);
    TAxis3D::GetPadAxis(gPad)->GetXaxis()->SetTitleSize(0.025);
    TAxis3D::GetPadAxis(gPad)->GetYaxis()->SetTitleSize(0.025);
    TAxis3D::GetPadAxis(gPad)->GetZaxis()->SetTitleSize(0.025);
    gPad->Modified();
    gPad->Update();

    std::cout<<"If the canvas has the focus, and if the mouse is on top of the canvas,"<<std::endl;
    std::cout<<"but not on any drawn object, the the following options are available:"<<std::endl;
    std::cout<<"X,Y,Z: move x,y, or z axis up by 25%"<<std::endl;
    std::cout<<"x,y,z: move x,y, or z axis down by 25%"<<std::endl;
    std::cout<<"P,T  : increase phi, or theta by 15 degrees"<<std::endl;
    std::cout<<"p,t  : decrease phi, or theta by 15 degrees"<<std::endl;
    std::cout<<"+,-  : zoom in or out"<<std::endl;
  }

  static void input()
  {
    int e=gPad->GetEvent(); //from http://root.cern.ch/root/html/tutorials/graphics/mandelbrot.C.html
    int key=gPad->GetEventX();
    if(e!=kKeyPress) return;
    switch(key)
    {
      case 'x':
      case 'y':
      case 'z': changerange(key-'x',true);
                break;
      case 'X':
      case 'Y':
      case 'Z': changerange(key-'X',false);
                break;
      case 'p':
      case 'P': changeangles(1,key=='p');
                break;
      case 't':
      case 'T': changeangles(2,key=='t');
                break;
      default : return;
    };
  }

  private:

  static void changerange(int axis, bool negativedirection)
  {
    double rmin[3],rmax[3];
    gPad->GetView()->GetRange(rmin,rmax);
    double difference=(rmax[axis]-rmin[axis])*0.25;
    if(negativedirection) difference*=-1;
    rmin[axis]+=difference;
    rmax[axis]+=difference;
    gPad->GetView()->SetRange(rmin,rmax);
    gPad->Modified();
    gPad->Update();
  }

  static void changeangles(int angle, bool negativedirection)
  {
    int irep;
    double phi=gPad->GetView()->GetLongitude();
    double theta=gPad->GetView()->GetLatitude();
    double difference=15;
    if(negativedirection) difference*=-1;
    if(angle==1) phi+=difference; else theta+=difference;
    gPad->GetView()->SetView(phi,theta,90,irep);
    gPad->SetPhi(-90-phi);
    gPad->SetTheta(90-theta);
    gPad->Modified();
    gPad->Update();
  }

  ClassDef(EventDisplayViewSetup,0);
};

}
#endif /* EventDisplay_src_dict_classes_EventDisplayViewSetup_h */
