/*
This script can be used to read root files which were produced by the event display.
Here is a sample file: http://home.fnal.gov/~ehrlich/events.root

You can run it like

root eventdisplay_reader.C\(\"events.root\"\)

or from inside of root with

.x eventdisplay_reader.C("events.root")
*/

TTree *tree;
TObjArray *objarray;
int nevents;
int currentevent;

void drawevent(int i)
{
  objarray->Clear();
  tree->GetEvent(i);
  objarray->Draw();
  gPad->Modified();
  gPad->Update();
}

void changerange(int axis, bool negativedirection)
{
  double rmin[3],rmax[3];
  gPad->GetView()->GetRange(rmin,rmax);
  double difference=(rmax[axis]-rmin[axis])*0.25;
  if(negativedirection) difference*=-1;
  rmin[axis]+=difference;
  rmax[axis]+=difference;
  gPad->GetView()->SetRange(rmin,rmax);
}

void changeangles(int angle, bool negativedirection)
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

void keyboard()
{
  int e=gPad->GetEvent(); //from http://root.cern.ch/root/html/tutorials/graphics/mandelbrot.C.html
  int key=gPad->GetEventX();
  if(e!=kKeyPress) return;
  switch(key)
  {
    case 'n': currentevent--;  //backward
              if(currentevent<0) currentevent=nevents-1;
              break;
    case 'm': currentevent++;  //forward
              if(currentevent>=nevents) currentevent=0;
              break;
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
  }
  drawevent(currentevent);
}

void setup()
{
  int irep;
  gPad->GetView()->SetView(180,0,90,irep);
  gPad->SetPhi(-90-180);
  gPad->SetTheta(90-0);
  gPad->GetView()->ShowAxis();
  gPad->SetFillColor(0);
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
  double rmin[3],rmax[3];
  gPad->GetView()->GetRange(rmin,rmax);
  rmin[0]=-800;
  rmax[0]=800;
  rmin[1]=-800;
  rmax[1]=800;
  gPad->GetView()->SetRange(rmin,rmax);
  gPad->Modified();
  gPad->Update();
  cout<<"If the canvas has the focus, and if the mouse is on top of the canvas,"<<endl; 
  cout<<"but not on any drawn object, the the following options are available:"<<endl;
  cout<<"m    : next event"<<endl;
  cout<<"n    : previous event"<<endl;
  cout<<"X,Y,Z: move x,y, or z axis up by 25%"<<endl;
  cout<<"x,y,z: move x,y, or z axis down by 25%"<<endl;
  cout<<"P,T  : increase phi, or theta by 15 degrees"<<endl;
  cout<<"p,t  : decrease phi, or theta by 15 degrees"<<endl;
  cout<<"+,-  : zoom in or out"<<endl;
}

void eventdisplay_reader(char* filename)
{
  TFile *file=new TFile(filename);
  tree=dynamic_cast<TTree*>(file->Get("Tree"));
  objarray=new TObjArray();
  objarray->SetOwner();
  tree->SetBranchAddress("Branch",&objarray);
  nevents=tree->GetEntries();
  if(nevents==0) return;
  currentevent=0;

  for(int i=0; i<20; i++)
  {
    float r,g,b;
    TColor *c;
    TColor::HLS2RGB(i*360/20,.5,.5,r,g,b);
    if(!gROOT->GetColor(i+2000)) c = new TColor(i+2000,r,g,b);
  }

  TCanvas *canvas=new TCanvas("Canvas","Mu2e Event",800,800);
  drawevent(0);
  setup();

//from http://root.cern.ch/root/html/tutorials/graphics/mandelbrot.C.html
  gPad->AddExec("KeyboardClick","keyboard()");
}
