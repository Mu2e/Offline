#include "DataInterface.h"
#include "Track.h"
#include "Straw.h"
#include "Cylinder.h"
#include "Cube.h"
#include "dict_classes/ComponentInfo.h"
#include "ContentSelector.h"

#include <TView.h>
#include <TAxis3D.h>
#include <TGFrame.h>
#include <TMath.h>

#include "FWCore/Framework/interface/Run.h"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "TargetGeom/inc/Target.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "HepPID/ParticleName.hh"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "ToyDP/inc/StrawHitCollection.hh"
#include "ToyDP/inc/SimParticleCollection.hh"
#include "ToyDP/inc/PointTrajectoryCollection.hh"
#include "ToyDP/inc/CaloCrystalHitCollection.hh"
#include "ToyDP/inc/CaloHitCollection.hh"
#include "ToyDP/inc/PhysicalVolumeInfoCollection.hh"

#include <TGeoVolume.h>

namespace mu2e_eventdisplay
{

DataInterface::DataInterface(const TGMainFrame *mainframe):
              _geometrymanager(NULL),_topvolume(NULL),_mainframe(mainframe),
              _showUnhitStraws(false), _showUnhitCrystals(false)
{
}

DataInterface::~DataInterface()
{
  removeAllComponents();
}

void DataInterface::startComponents()
{
  std::list<boost::shared_ptr<VirtualShape> >::const_iterator iter;
  for(iter=_components.begin(); iter!=_components.end(); iter++)
  {
    (*iter)->start();
  }
}

void DataInterface::updateComponents(double time)
{
  std::vector<boost::shared_ptr<Track> >::const_iterator track;
  for(track=_tracks.begin(); track!=_tracks.end(); track++)
  {
    (*track)->update(time);
  }
  if(_showUnhitStraws)
  {
    std::map<int,boost::shared_ptr<Straw> >::const_iterator straw;
    for(straw=_straws.begin(); straw!=_straws.end(); straw++)
    {
      straw->second->update(time);
    }
  }
  else
  {
    std::vector<boost::shared_ptr<Straw> >::const_iterator hit;
    for(hit=_hits.begin(); hit!=_hits.end(); hit++)
    {
      (*hit)->update(time);
    }
  }
  if(_showUnhitCrystals)
  {
    std::map<int,boost::shared_ptr<Cube> >::const_iterator crystal;
    for(crystal=_crystals.begin(); crystal!=_crystals.end(); crystal++)
    {
      crystal->second->update(time);
    }
  }
  else
  {
    std::vector<boost::shared_ptr<Cube> >::const_iterator crystalhit;
    for(crystalhit=_crystalhits.begin(); crystalhit!=_crystalhits.end(); crystalhit++)
    {
      (*crystalhit)->update(time);
    }
  }
}

void DataInterface::createGeometryManager()
{
  _geometrymanager = new TGeoManager("GeoManager", "GeoManager");
//  _geometrymanager->SetVerboseLevel(0);     //wait till root version 5.22 or so
  _topvolume = _geometrymanager->MakeBox("TopVolume", NULL, 1000, 1000, 1500);
  _geometrymanager->SetTopVolume(_topvolume);
  _geometrymanager->SetTopVisible(false);
  _geometrymanager->CloseGeometry();
  _geometrymanager->SetVisLevel(4);
  _topvolume->SetVisibility(0);
  _topvolume->SetLineColor(0);
  _topvolume->Draw("ogle");
  int irep=0;
  gPad->GetView()->SetView(180,65,90,irep);
  gPad->SetPhi(-90-180);
  gPad->SetTheta(90-65);
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
}

void DataInterface::fillGeometry()
{
  removeAllComponents();
  createGeometryManager();
  resetBoundaryP(_trackerMinmax);
  resetBoundaryP(_targetMinmax);
  resetBoundaryP(_calorimeterMinmax);
  resetBoundaryP(_tracksMinmax);

  edm::Service<mu2e::GeometryService> geom;

  const mu2e::SimpleConfig &config = geom->config();
  _xOffset=config.getDouble("mu2e.solenoidOffset");    //between Mu2e and Tracker coordinates
  _zOffset=-config.getDouble("mu2e.detectorSystemZ0"); //between Mu2e and Tracker coordinates
  _zOffsetDS=1800.0;                                   //between DS and Tracker coordinates
//coordinate transformation between "World" and Mu2e coordinates
//(taken from CosmicRaysShieldMaker class)
  std::vector<double> worldHLen;
  config.getVectorDouble("world.halfLengths", worldHLen, 3);
  double floorThick =  config.getDouble("hall.floorThick");
  double yFloor     = -worldHLen[1] + floorThick;
  CLHEP::Hep3Vector _mu2eOriginInWorld = CLHEP::Hep3Vector(config.getDouble("world.mu2eOrigin.xoffset"),
                                                           config.getDouble("world.mu2eOrigin.height") + yFloor,
                                                           config.getDouble("world.mu2eOrigin.zoffset"));

  if(geom->hasElement<mu2e::TTracker>())
  {
//Straws
    mu2e::GeomHandle<mu2e::TTracker> ttracker;  
    const std::deque<mu2e::Straw>& allStraws = ttracker->getAllStraws();
    std::deque<mu2e::Straw>::const_iterator iter;
    for(iter=allStraws.begin(); iter!=allStraws.end(); iter++)
    {
      const mu2e::Straw &s=*iter;
      const CLHEP::Hep3Vector& p = s.getMidPoint();
      const CLHEP::Hep3Vector& d = s.getDirection();
      double x = p.x();
      double y = p.y();
      double z = p.z();
      double theta = d.theta();
      double phi = d.phi();
//      double r = s.getRadius();
      double l = s.getHalfLength();
      int idStraw =  s.Id().getStraw();
      int idLayer =  s.Id().getLayer();
      int idSector =  s.Id().getSector();
      int idDevice =  s.Id().getDevice();
      int index = s.Index().asInt();

      char c[200];
      sprintf(c,"Straw %i  Layer %i  Sector %i  Device %i",idStraw,idLayer,idSector,idDevice);
      boost::shared_ptr<ComponentInfo> info(new ComponentInfo());
      info->setName(c);
      info->setText(0,c);
      boost::shared_ptr<Straw> shape(new Straw(x,y,z, NAN, theta, phi, l, _geometrymanager, _topvolume, _mainframe, info, false));
      _components.push_back(shape);
      _straws[index]=shape;
    }

//Support Structure
    double innerRadius=ttracker->getSupportParams().innerRadius;
    double outerRadius=ttracker->getSupportParams().outerRadius;
    double zHalfLength=ttracker->getTrackerEnvelopeParams().zHalfLength;
    findBoundaryP(_trackerMinmax, outerRadius, outerRadius, zHalfLength);
    findBoundaryP(_trackerMinmax, -outerRadius, -outerRadius, -zHalfLength);

    char c[200];
    boost::shared_ptr<ComponentInfo> info(new ComponentInfo());
    sprintf(c,"TTracker Support Structure");
    info->setName(c);
    info->setText(0,c);
    sprintf(c,"Inner Radius %.f mm  Outer Radius %.f mm",innerRadius/CLHEP::mm,outerRadius/CLHEP::mm);
    info->setText(1,c);
    sprintf(c,"Length %.f mm",2.0*zHalfLength/CLHEP::mm);
    info->setText(2,c);
    sprintf(c,"Center at x: 0 mm, y: 0 mm, z: 0 mm");
    info->setText(3,c);
    boost::shared_ptr<Cylinder> shape(new Cylinder(0,0,0, 0,0,0, 
                                          zHalfLength,innerRadius,outerRadius,
                                          _geometrymanager, _topvolume, _mainframe, info, true));
    _components.push_back(shape);
    _supportstructures.push_back(shape);

//Envelope
    innerRadius=ttracker->getTrackerEnvelopeParams().innerRadius;
    outerRadius=ttracker->getTrackerEnvelopeParams().outerRadius;
    zHalfLength=ttracker->getTrackerEnvelopeParams().zHalfLength;

    boost::shared_ptr<ComponentInfo> infoEnvelope(new ComponentInfo());
    sprintf(c,"TTracker Envelope");
    infoEnvelope->setName(c);
    infoEnvelope->setText(0,c);
    sprintf(c,"Inner Radius %.f mm  Outer Radius %.f mm",innerRadius/CLHEP::mm,outerRadius/CLHEP::mm);
    infoEnvelope->setText(1,c);
    sprintf(c,"Length %.f mm",2.0*zHalfLength/CLHEP::mm);
    infoEnvelope->setText(2,c);
    sprintf(c,"Center at x: 0 mm, y: 0 mm, z: 0 mm");
    infoEnvelope->setText(3,c);
    boost::shared_ptr<Cylinder> shapeEnvelope(new Cylinder(0,0,0, 0,0,0, 
                                                  zHalfLength,innerRadius,outerRadius,
                                                  _geometrymanager, _topvolume, _mainframe, infoEnvelope, true));
    _components.push_back(shapeEnvelope);
    _supportstructures.push_back(shapeEnvelope);

//ToyDS
    innerRadius=config.getDouble("toyDS.rIn");
    outerRadius=config.getDouble("toyDS.rOut");
    zHalfLength=config.getDouble("toyDS.halfLength");
    double z=config.getDouble("toyDS.z0")+_zOffset;

    boost::shared_ptr<ComponentInfo> infoToyDS(new ComponentInfo());
    sprintf(c,"Toy DS");
    infoToyDS->setName(c);
    infoToyDS->setText(0,c);
    sprintf(c,"Inner Radius %.f mm  Outer Radius %.f mm",innerRadius/CLHEP::mm,outerRadius/CLHEP::mm);
    infoToyDS->setText(1,c);
    sprintf(c,"Length %.f mm",2.0*zHalfLength/CLHEP::mm);
    infoToyDS->setText(2,c);
    sprintf(c,"Center at x: 0 mm, y: 0 mm, z: %.f mm",z/CLHEP::mm);
    infoToyDS->setText(3,c);
    boost::shared_ptr<Cylinder> shapeToyDS(new Cylinder(0,0,z, 0,0,0, 
                                               zHalfLength,innerRadius,outerRadius,
                                               _geometrymanager, _topvolume, _mainframe, infoToyDS, true));
    _components.push_back(shapeToyDS);
    _otherstructures.push_back(shapeToyDS);
  }

  if(geom->hasElement<mu2e::Target>())
  {
    mu2e::GeomHandle<mu2e::Target> target;  
    double radius=target->cylinderRadius();
    double length=target->cylinderLength();
    double z=target->cylinderCenter()+_zOffsetDS;
    unsigned int n=target->nFoils();
    for(unsigned int i=0; i<n; i++)
    {
//      const TargetFoil &foil=target->foil(i);
    }

    findBoundaryP(_targetMinmax, radius, radius, z+length);
    findBoundaryP(_targetMinmax, -radius, -radius, z-length);

    char c[200];
    boost::shared_ptr<ComponentInfo> info(new ComponentInfo());
    sprintf(c,"Target");
    info->setName(c);
    info->setText(0,c);
    sprintf(c,"Outer Radius %.f mm",radius/CLHEP::mm);
    info->setText(1,c);
    sprintf(c,"Length %.f mm",length/CLHEP::mm);
    info->setText(2,c);
    sprintf(c,"Center at x: 0 mm, y: 0 mm, z: %.f mm",z/CLHEP::mm);
    info->setText(3,c);
    boost::shared_ptr<Cylinder> shape(new Cylinder(0,0,z, 0,0,0, length/2.0,0,radius,
                                          _geometrymanager, _topvolume, _mainframe, info, true));
    _components.push_back(shape);
    _supportstructures.push_back(shape);
  }
  
  if(geom->hasElement<mu2e::Calorimeter>())
  {
    mu2e::GeomHandle<mu2e::Calorimeter> calo;  
    unsigned int n=calo->nVane();
    for(unsigned int i=0; i<n; i++)
    { 
      const mu2e::Vane &v=calo->getVane(i);
      double x=v.getOrigin().x()+_xOffset;
      double y=v.getOrigin().y();
      double z=v.getOrigin().z()+_zOffset;
      int    id=v.Id();
      double sx=v.getSize().x();
      double sy=v.getSize().y();
      double sz=v.getSize().z();
      double phi=v.getRotation()->phi();
      double theta=v.getRotation()->theta();
      double psi=v.getRotation()->psi();

      findBoundaryP(_calorimeterMinmax, x+sx, y+sy, z+sz);
      findBoundaryP(_calorimeterMinmax, x-sx, y-sy, z-sz);

      char c[200];
      boost::shared_ptr<ComponentInfo> info(new ComponentInfo());
      sprintf(c,"Vane %i",id);
      info->setName(c);
      info->setText(0,c);
      sprintf(c,"Dimension  ∆x: %.f mm, ∆y: %.f mm, ∆z: %.f mm",sx/CLHEP::mm,sy/CLHEP::mm,sz/CLHEP::mm);
      info->setText(1,c);
      sprintf(c,"Rotation phi: %.f °, theta: %.f °, psi: %.f °",phi/CLHEP::deg,theta/CLHEP::deg,psi/CLHEP::deg);
      info->setText(2,c);
      sprintf(c,"Center at x: %.f mm, y: %.f mm, z: %.f mm",x/CLHEP::mm,y/CLHEP::mm,z/CLHEP::mm);
      info->setText(3,c);
      boost::shared_ptr<Cube> shape(new Cube(x,y,z,  sx,sy,sz,  phi,theta,psi,   NAN,0, 
                                        _geometrymanager, _topvolume, _mainframe, info, true));
      _components.push_back(shape);
      _supportstructures.push_back(shape);
    }

    unsigned int roPerCrystal=calo->nROPerCrystal();
    unsigned int nro=calo->nRO();
    for(unsigned int i=0; i<nro; i+=roPerCrystal)
    {
      int vaneid=calo->getVaneByRO(i);
      int crystalid=calo->getCrystalByRO(i);
      int rPos=calo->getCrystalRByRO(i);
      int zPos=calo->getCrystalZByRO(i);
      double crystalHalfSize=calo->crystalHalfSize();

      const mu2e::Vane &v=calo->getVane(vaneid);
      double x=v.getOrigin().x()+_xOffset;
      double y=v.getOrigin().y();
      double z=v.getOrigin().z()+_zOffset;
      double theta=v.getRotation()->theta();
      double phi=v.getRotation()->phi();
      double psi=v.getRotation()->psi();
      double sx=v.getSize().x();
      double sy=v.getSize().y();
      double sz=v.getSize().z();

      //Start with an unrotated vane centered at (0,0,0).
      //Before the rotation, the vector from the center of the vane 
      //to the center of a crystal is (crystalX,crystalY,crystalZ).
      double crystalX=0;
      double crystalY=-sy+crystalHalfSize*(2.0*rPos+1.0);
      double crystalZ=-sz+crystalHalfSize*(2.0*zPos+1.0);

      double st=sin(theta);
      double ct=cos(theta);
      double sp=sin(phi);
      double cp=cos(phi);
      double ss=sin(psi);
      double cs=cos(psi);

      //After the rotation of the vane, the vector from the center of the vane 
      //to the center of a crystal is (rotatedX,rotatedY,rotatedZ).
      double rotatedX,rotatedY,rotatedZ;
      VirtualShape::rotate(crystalX,crystalY,crystalZ,  rotatedX,rotatedY,rotatedZ,  sp,cp,st,ct,ss,cs);

      //After the vane gets shifted from (0,0,0) to (x,y,z), 
      //the crystal centers need to be shifted to (x+rotatedX,y+rotatedY,z+rotatedZ).
      char c[200];
      boost::shared_ptr<ComponentInfo> info(new ComponentInfo());
      sprintf(c,"Vane %i, Crystal %i",vaneid,crystalid);
      info->setName(c);
      info->setText(0,c);
      sprintf(c,"Dimension  ∆x: %.f mm, ∆y: %.f mm, ∆z: %.f mm",sx/CLHEP::mm,crystalHalfSize/CLHEP::mm,crystalHalfSize/CLHEP::mm);
      info->setText(1,c);
      sprintf(c,"Rotation phi: %.f °, theta: %.f °, psi: %.f °",phi/CLHEP::deg,theta/CLHEP::deg,psi/CLHEP::deg);
      info->setText(2,c);
      sprintf(c,"Center at x: %.f mm, y: %.f mm, z: %.f mm",(x+rotatedX)/CLHEP::mm,(y+rotatedY)/CLHEP::mm,(z+rotatedZ)/CLHEP::mm);
      info->setText(3,c);
      boost::shared_ptr<Cube> shape(new Cube(x+rotatedX,y+rotatedY,z+rotatedZ,  sx,crystalHalfSize,crystalHalfSize,  
                                        phi,theta,psi,  NAN,0,  
                                        _geometrymanager, _topvolume, _mainframe, info, false));
      _components.push_back(shape);
      _crystals[crystalid]=shape;
    }
  }

  if(geom->hasElement<mu2e::CosmicRayShield>())
  {
    mu2e::GeomHandle<mu2e::CosmicRayShield> crs;  
    std::string steelshieldnames[6]={"CRSSteelTopShield",
                                     "CRSSteelBottomShield",
                                     "CRSSteelLeftShield",
                                     "CRSSteelRightShield",
                                     "CRSSteelBackShield",
                                     "CRSSteelFrontShield"};
    for(int i=0; i<6; i++)
    {
//TODO: CosmicRayShield needs to have a public function 
//which returns the entire map so that one can interate over it
//without knowing the names of the elements
      const mu2e::CosmicRayShieldSteelShield &steelshield=crs->getCosmicRayShieldSteelShield(steelshieldnames[i].c_str());
      double x=steelshield.getGlobalOffset().x()-_mu2eOriginInWorld.x()+_xOffset;
      double y=steelshield.getGlobalOffset().y()-_mu2eOriginInWorld.y();
      double z=steelshield.getGlobalOffset().z()-_mu2eOriginInWorld.z()+_zOffset;
      double dx=steelshield.getHalfLengths()[0];
      double dy=steelshield.getHalfLengths()[1];
      double dz=steelshield.getHalfLengths()[2];
      char c[200];
      boost::shared_ptr<ComponentInfo> info(new ComponentInfo());
      sprintf(c,"Cosmic Ray Steel Shield %i (%s)",i,steelshieldnames[i].c_str());
      info->setName(c);
      info->setText(0,c);
      sprintf(c,"Dimension  ∆x: %.f mm, ∆y: %.f mm, ∆z: %.f mm",dx/CLHEP::mm,dy/CLHEP::mm,dz/CLHEP::mm);
      info->setText(1,c);
      sprintf(c,"Center at x: %.f mm, y: %.f mm, z: %.f mm",x/CLHEP::mm,y/CLHEP::mm,z/CLHEP::mm);
      info->setText(2,c);

      double holeRadius=steelshield.getHoleRadius();
      if(holeRadius==0)
      {
        boost::shared_ptr<Cube> shape(new Cube(x,y,z,  dx,dy,dz,  0,0,0, NAN,0,
                                          _geometrymanager, _topvolume, _mainframe, info, true));
        _components.push_back(shape);
        _otherstructures.push_back(shape);
      }
      else
      {
//TODO: This needs to be replaced by a shape with a hole
        boost::shared_ptr<Cube> shape(new Cube(x,y,z,  dx,dy,dz,  0,0,0, NAN,0,
                                          _geometrymanager, _topvolume, _mainframe, info, true));
        _components.push_back(shape);
        _otherstructures.push_back(shape);
      }
    }
  }
}

void DataInterface::makeSupportStructuresVisible(bool visible)
{
  std::vector<boost::shared_ptr<VirtualShape> >::const_iterator structure;
  for(structure=_supportstructures.begin(); structure!=_supportstructures.end(); structure++)
  {
    (*structure)->setDefaultVisibility(visible);
    (*structure)->start();
  }

  //tracks and straws don't have to be pushed into the foreground if the structure is removed
  if(visible) toForeground();
} 

void DataInterface::makeOtherStructuresVisible(bool visible)
{
  std::vector<boost::shared_ptr<VirtualShape> >::const_iterator structure;
  for(structure=_otherstructures.begin(); structure!=_otherstructures.end(); structure++)
  {
    (*structure)->setDefaultVisibility(visible);
    (*structure)->start();
  }

  //tracks and straws don't have to be pushed into the foreground if the structure is removed
  if(visible) toForeground();
} 

void DataInterface::makeStrawsVisibleBeforeStart(bool visible)
{
  std::map<int,boost::shared_ptr<Straw> >::const_iterator straw;
  for(straw=_straws.begin(); straw!=_straws.end(); straw++)
  {
    straw->second->setDefaultVisibility(visible);
    straw->second->start();
  }
  _showUnhitStraws=visible;
} 

void DataInterface::makeCrystalsVisibleBeforeStart(bool visible)
{
  std::map<int,boost::shared_ptr<Cube> >::const_iterator crystal;
  for(crystal=_crystals.begin(); crystal!=_crystals.end(); crystal++)
  {
    crystal->second->setDefaultVisibility(visible);
    crystal->second->start();
  }
  _showUnhitCrystals=visible;
} 
  
void DataInterface::toForeground()
{
  std::map<int,boost::shared_ptr<Straw> >::const_iterator straw;
  for(straw=_straws.begin(); straw!=_straws.end(); straw++)
  {
    straw->second->toForeground();
  }

  std::map<int,boost::shared_ptr<Cube> >::const_iterator crystal;
  for(crystal=_crystals.begin(); crystal!=_crystals.end(); crystal++)
  {
    crystal->second->toForeground();
  }

  std::vector<boost::shared_ptr<Track> >::const_iterator track;
  for(track=_tracks.begin(); track!=_tracks.end(); track++)
  {
    (*track)->toForeground();
  }
}

void DataInterface::useHitColors(bool hitcolors, bool whitebackground)
{
  std::vector<boost::shared_ptr<Straw> >::const_iterator hit;
  for(hit=_hits.begin(); hit!=_hits.end(); hit++)
  {
    double time=(*hit)->getStartTime();
    if(hitcolors)
    {
      int color=TMath::FloorNint(20.0*(time-_hitsTimeMinmax.mint)/(_hitsTimeMinmax.maxt-_hitsTimeMinmax.mint));
      if(color>=20) color=19;
      if(color<=0) color=0;
      color+=2000;
      (*hit)->setColor(color);
    }
    else (*hit)->setColor(whitebackground?1:0);
  }
  std::vector<boost::shared_ptr<Cube> >::const_iterator crystalhit;
  for(crystalhit=_crystalhits.begin(); crystalhit!=_crystalhits.end(); crystalhit++)
  {
    double time=(*crystalhit)->getStartTime();
    if(hitcolors)
    {
      int color=TMath::FloorNint(20.0*(time-_hitsTimeMinmax.mint)/(_hitsTimeMinmax.maxt-_hitsTimeMinmax.mint));
      if(color>=20) color=19;
      if(color<=0) color=0;
      color+=2000;
      (*crystalhit)->setColor(color);
    }
    else (*crystalhit)->setColor(whitebackground?1:0);
  }
}

void DataInterface::useTrackColors(bool trackcolors, bool whitebackground)
{
  std::vector<boost::shared_ptr<Track> >::const_iterator track;
  for(track=_tracks.begin(); track!=_tracks.end(); track++)
  {
    if(trackcolors)
    {
      int particleid=(*track)->getParticleId();
      int color=kGray;
      switch(particleid)
      {
        case   11:
        case  -11: color=2; break;   //e+,e-
        case   13:
        case  -13: color=3; break;   //mu+,mu-
        case   22: color=4; break;   //gamma
        case 2112: color=6; break;   //n0
        case   12: 
        case  -12: 
        case   14: 
        case  -14: 
        case   16: 
        case  -16: color=28; break;   //neutrinos
      };
      (*track)->setColor(color);
    }
    else (*track)->setColor(whitebackground?1:0);
  }
}

void DataInterface::resetBoundaryT(timeminmax &m)
{
  m.mint=NAN;
  m.maxt=NAN;
}

void DataInterface::resetBoundaryP(spaceminmax &m)
{
  m.minx=NAN;
  m.miny=NAN;
  m.minz=NAN;
  m.maxx=NAN;
  m.maxy=NAN;
  m.maxz=NAN;
}

DataInterface::spaceminmax DataInterface::getSpaceBoundary(bool useTarget, bool useCalorimeter, bool useTracks)
{
  spaceminmax m;
  resetBoundaryP(m);
  findBoundaryP(m, _trackerMinmax.minx, _trackerMinmax.miny, _trackerMinmax.minz);
  findBoundaryP(m, _trackerMinmax.maxx, _trackerMinmax.maxy, _trackerMinmax.maxz);
  if(useTarget)
  {
    findBoundaryP(m, _targetMinmax.minx, _targetMinmax.miny, _targetMinmax.minz);
    findBoundaryP(m, _targetMinmax.maxx, _targetMinmax.maxy, _targetMinmax.maxz);
  }
  if(useCalorimeter)
  {
    findBoundaryP(m, _calorimeterMinmax.minx, _calorimeterMinmax.miny, _calorimeterMinmax.minz);
    findBoundaryP(m, _calorimeterMinmax.maxx, _calorimeterMinmax.maxy, _calorimeterMinmax.maxz);
  }
  if(useTracks)
  {
    findBoundaryP(m, _tracksMinmax.minx, _tracksMinmax.miny, _tracksMinmax.minz);
    findBoundaryP(m, _tracksMinmax.maxx, _tracksMinmax.maxy, _tracksMinmax.maxz);
  }

  if(isnan(m.minx)) m.minx=-1000;
  if(isnan(m.miny)) m.miny=-1000;
  if(isnan(m.minz)) m.minz=-1000;
  if(isnan(m.maxx)) m.maxx=1000;
  if(isnan(m.maxy)) m.maxy=1000;
  if(isnan(m.maxz)) m.maxz=1000;
  return m;
}

void DataInterface::findBoundaryT(timeminmax &m, double t)
{
  if(isnan(m.mint) || t<m.mint) m.mint=t;
  if(isnan(m.maxt) || t>m.maxt) m.maxt=t;
}

void DataInterface::findBoundaryP(spaceminmax &m, double x, double y, double z)
{
  if(isnan(m.minx) || x<m.minx) m.minx=x;
  if(isnan(m.miny) || y<m.miny) m.miny=y;
  if(isnan(m.minz) || z<m.minz) m.minz=z;
  if(isnan(m.maxx) || x>m.maxx) m.maxx=x;
  if(isnan(m.maxy) || y>m.maxy) m.maxy=y;
  if(isnan(m.maxz) || z>m.maxz) m.maxz=z;
}

void DataInterface::fillEvent(const ContentSelector *contentSelector)
{
  removeNonGeometryComponents();
  if(!_geometrymanager) createGeometryManager();
  resetBoundaryT(_hitsTimeMinmax);
  resetBoundaryT(_tracksTimeMinmax);

  _numberHits=0;
  _numberCrystalHits=0;

  const mu2e::StepPointMCCollection *steppointMChits=contentSelector->getSelectedHitCollection<mu2e::StepPointMCCollection>();
  if(steppointMChits!=NULL)  
  {
    _numberHits=steppointMChits->size();
    std::vector<mu2e::StepPointMC>::const_iterator iter;
    for(iter=steppointMChits->begin(); iter!=steppointMChits->end(); iter++)
    {
      const mu2e::StepPointMC& hit = *iter;
      int strawindex = hit.strawIndex().asInt();
      int trackid = hit.trackId().asInt();
      double time = hit.time();
      double energy = hit.eDep();
      std::map<int,boost::shared_ptr<Straw> >::iterator straw=_straws.find(strawindex);
      if(straw!=_straws.end() && !isnan(time))
      {
        double previousStartTime=straw->second->getStartTime();
        if(isnan(previousStartTime))
        {
          findBoundaryT(_hitsTimeMinmax, time);  //is it Ok to exclude all following hits from the time window?
          straw->second->setStartTime(time);
          straw->second->start();
          char c[100];
          sprintf(c,"hit time(s): %gns",time/CLHEP::ns);
          straw->second->getComponentInfo()->setText(1,c);
          sprintf(c,"deposited energy(s): %geV",energy/CLHEP::eV);
          straw->second->getComponentInfo()->setText(2,c);
          sprintf(c,"track id(s): %i",trackid);
          straw->second->getComponentInfo()->setText(3,c);
          _hits.push_back(straw->second);
        }
        else
        {
          straw->second->getComponentInfo()->expandLine(1,time/CLHEP::ns,"ns");
          straw->second->getComponentInfo()->expandLine(2,energy/CLHEP::eV,"eV");
          straw->second->getComponentInfo()->expandLine(3,trackid,"");
        }
      }
    }
  }

  const mu2e::StrawHitCollection *strawhits=contentSelector->getSelectedHitCollection<mu2e::StrawHitCollection>();
  if(strawhits!=NULL)  
  {
    _numberHits=strawhits->size();
    std::vector<mu2e::StrawHit>::const_iterator iter;
    for(iter=strawhits->begin(); iter!=strawhits->end(); iter++)
    {
      const mu2e::StrawHit& hit = *iter;
      int strawindex = hit.strawIndex().asInt();
      double time = hit.time();
      double dt = hit.dt();
      double energy = hit.energyDep();
      std::map<int,boost::shared_ptr<Straw> >::iterator straw=_straws.find(strawindex);
      if(straw!=_straws.end() && !isnan(time))
      {
        double previousStartTime=straw->second->getStartTime();
        if(isnan(previousStartTime))
        {
          findBoundaryT(_hitsTimeMinmax, time);  //is it Ok to exclude all following hits from the time window?
          straw->second->setStartTime(time);
          straw->second->start();
          char c[100];
          sprintf(c,"hit time(s): %gns",time/CLHEP::ns);
          straw->second->getComponentInfo()->setText(1,c);
          sprintf(c,"deposited energy(s): %geV",energy/CLHEP::eV);
          straw->second->getComponentInfo()->setText(2,c);
          sprintf(c,"hit time interval(s): %gns",dt/CLHEP::ns);
          straw->second->getComponentInfo()->setText(3,c);
          _hits.push_back(straw->second);
        }
        else
        {
          straw->second->getComponentInfo()->expandLine(1,time/CLHEP::ns,"ns");
          straw->second->getComponentInfo()->expandLine(2,energy/CLHEP::eV,"eV");
          straw->second->getComponentInfo()->expandLine(3,dt/CLHEP::ns,"ns");
        }
      }
    }
  }


  const mu2e::CaloCrystalHitCollection *calocrystalhits=contentSelector->getSelectedCaloHitCollection<mu2e::CaloCrystalHitCollection>();
  if(calocrystalhits!=NULL)  
  {
    _numberCrystalHits=calocrystalhits->size();
    std::vector<mu2e::CaloCrystalHit>::const_iterator iter;
    for(iter=calocrystalhits->begin(); iter!=calocrystalhits->end(); iter++)
    {
      const mu2e::CaloCrystalHit& calohit = *iter;
      int crystalid = calohit.id();
      double time = calohit.time();
      double energy = calohit.energyDep();
      std::map<int,boost::shared_ptr<Cube> >::iterator crystal=_crystals.find(crystalid);
      if(crystal!=_crystals.end() && !isnan(time))
      {
        double previousStartTime=crystal->second->getStartTime();
        if(isnan(previousStartTime))
        {
          findBoundaryT(_hitsTimeMinmax, time);  //is it Ok to exclude all following hits from the time window?
          crystal->second->setStartTime(time);
          crystal->second->start();
          char c[100];
          sprintf(c,"hit time(s): %gns",time/CLHEP::ns);
          crystal->second->getComponentInfo()->setText(1,c);
          sprintf(c,"deposited energy(s): %geV",energy/CLHEP::eV);
          crystal->second->getComponentInfo()->setText(2,c);
          _crystalhits.push_back(crystal->second);
        }
        else
        {
          crystal->second->getComponentInfo()->expandLine(1,time/CLHEP::ns,"ns");
          crystal->second->getComponentInfo()->expandLine(2,energy/CLHEP::eV,"eV");
        }
      }
    }
  }

  const mu2e::CaloHitCollection *calohits=contentSelector->getSelectedCaloHitCollection<mu2e::CaloHitCollection>();
  edm::Service<mu2e::GeometryService> geoservice;
  if(calohits!=NULL && geoservice->hasElement<mu2e::Calorimeter>())
  {
    mu2e::GeomHandle<mu2e::Calorimeter> calo;  
    _numberCrystalHits=calohits->size();  //this is not accurate since the return value gives the RO hits
    std::vector<mu2e::CaloHit>::const_iterator iter;
    for(iter=calohits->begin(); iter!=calohits->end(); iter++)
    {
      const mu2e::CaloHit& calohit = *iter;
      int roid = calohit.id();
      int crystalid=calo->getCrystalByRO(roid);
      double time = calohit.time();
      double energy = calohit.energyDep();
      std::map<int,boost::shared_ptr<Cube> >::iterator crystal=_crystals.find(crystalid);
      if(crystal!=_crystals.end() && !isnan(time))
      {
        double previousStartTime=crystal->second->getStartTime();
        if(isnan(previousStartTime))
        {
          findBoundaryT(_hitsTimeMinmax, time);  //is it Ok to exclude all following hits from the time window?
          crystal->second->setStartTime(time);
          crystal->second->start();
          char c[100];
          sprintf(c,"hit time(s): %gns",time/CLHEP::ns);
          crystal->second->getComponentInfo()->setText(1,c);
          sprintf(c,"deposited energy(s): %geV",energy/CLHEP::eV);
          crystal->second->getComponentInfo()->setText(2,c);
          sprintf(c,"RO ID(s): %i",roid);
          crystal->second->getComponentInfo()->setText(3,c);
          _crystalhits.push_back(crystal->second);
        }
        else
        {
          crystal->second->getComponentInfo()->expandLine(1,time/CLHEP::ns,"ns");
          crystal->second->getComponentInfo()->expandLine(2,energy/CLHEP::eV,"eV");
          crystal->second->getComponentInfo()->expandLine(3,roid,"");
        }
      }
    }
  }

  unsigned int physicalVolumeEntries=0;
  const mu2e::PhysicalVolumeInfoCollection *physicalVolumes=contentSelector->getPhysicalVolumeInfoCollection();
  if(physicalVolumes!=NULL)
  {
    physicalVolumeEntries=physicalVolumes->size();
  }

  resetBoundaryP(_tracksMinmax);
  std::vector<const mu2e::SimParticleCollection*> trackCollectionVector=contentSelector->getSelectedTrackCollection<mu2e::SimParticleCollection>();
  for(unsigned int i=0; i<trackCollectionVector.size(); i++) 
  {
    const mu2e::SimParticleCollection *simParticles=trackCollectionVector[i];
    MapVector<mu2e::SimParticle>::const_iterator iter;
    for(iter=simParticles->begin(); iter!=simParticles->end(); iter++)
    {
      const mu2e::SimParticle& particle = iter->second;
      int id = particle.id().asInt();
      int parentid = particle.parentId().asInt();
      int particleid=particle.pdgId();
      std::string particlename=HepPID::particleName(particle.pdgId());
      unsigned int startVolume=particle.startVolumeIndex();
      unsigned int endVolume  =particle.endVolumeIndex();
      std::string startVolumeName="unknown volume";
      std::string endVolumeName="unknown volume";
      if(startVolume<physicalVolumeEntries && startVolume>=0) startVolumeName=physicalVolumes->at(startVolume).name();     
      if(endVolume<physicalVolumeEntries && endVolume>=0) endVolumeName=physicalVolumes->at(endVolume).name();     

      double x1=particle.startPosition().x()+_xOffset;
      double y1=particle.startPosition().y();
      double z1=particle.startPosition().z()+_zOffset;
      double t1=particle.startGlobalTime();
      double e1=particle.startMomentum().e();
      double x2=particle.endPosition().x()+_xOffset;
      double y2=particle.endPosition().y();
      double z2=particle.endPosition().z()+_zOffset;
      double t2=particle.endGlobalTime();
      double e2=particle.endMomentum().e();
      findBoundaryT(_tracksTimeMinmax, t1);
      findBoundaryT(_tracksTimeMinmax, t2);
      findBoundaryP(_tracksMinmax, x1, y1, z1);
      findBoundaryP(_tracksMinmax, x2, y2, z2);

      char c0[200], c1[200],c2[200],c3[200],c4[200];
      sprintf(c0,"Track %i  %s",id,particlename.c_str());
      if(parentid!=0) sprintf(c1,"Track %i  %s  Parent %i",id,particlename.c_str(),parentid);
      else sprintf(c1,"Track %i  %s  generated track",id,particlename.c_str());
      sprintf(c2,"Start Energy %gMeV  End Energy %gMeV",e1/CLHEP::MeV,e2/CLHEP::MeV);
      sprintf(c3,"Created by %s in %s",particle.creationCode().name().c_str(),startVolumeName.c_str());
      sprintf(c4,"Destroyed by %s in %s",particle.stoppingCode().name().c_str(),endVolumeName.c_str());
      boost::shared_ptr<ComponentInfo> info(new ComponentInfo());
      info->setName(c0);
      info->setText(0,c1);
      info->setText(1,c2);
      info->setText(2,c3);
      info->setText(3,c4);
      info->setText(4,"Daughter IDs:");
      std::vector<MapVectorKey>::const_iterator daughter;
      for(daughter=particle.daughterIds().begin(); 
          daughter!=particle.daughterIds().end();
          daughter++)
      {
        info->expandLine(4,daughter->asInt(),"");
      }
      boost::shared_ptr<Track> shape(new Track(x1,y1,z1,t1, x2,y2,z2,t2, particleid, _geometrymanager, _topvolume, _mainframe, info));
      findTrajectory(contentSelector,shape,id);
      _components.push_back(shape);
      _tracks.push_back(shape);
    }
  }
}

bool DataInterface::findTrajectory(const ContentSelector *contentSelector,
                                   boost::shared_ptr<Track> track, int id)
{
  const mu2e::PointTrajectoryCollection *pointTrajectories=contentSelector->getPointTrajectoryCollection();
  if(pointTrajectories!=NULL)
  {
    MapVector<mu2e::PointTrajectory>::const_iterator iter;
    for(iter=pointTrajectories->begin(); iter!=pointTrajectories->end(); iter++)
    {
      const mu2e::PointTrajectory& trajectory = iter->second;
      int trajectory_id = trajectory.simId();
      if(id==trajectory_id)
      {
        const std::vector<CLHEP::Hep3Vector>& p_vec=trajectory.points();
        for(unsigned int i=0; i<p_vec.size(); i++)
        {
//TODO: no time information available from PointTrajectory
//the assumption that the time period between each 
//recorded trajectory point is equal is not valid
          track->addTrajectoryPoint(p_vec[i].x(),p_vec[i].y(),p_vec[i].z()+_zOffsetDS,i,p_vec.size());
        }
        return true;
      }
    }
  }
  return false;
}

void DataInterface::removeNonGeometryComponents()
{
  std::list<boost::shared_ptr<VirtualShape> >::iterator iter=_components.begin();
  while(iter!=_components.end()) 
  {
    if(!(*iter)->isGeometry()) {iter=_components.erase(iter);}
    else 
    {
      (*iter)->setStartTime(NAN);
      (*iter)->start(); 
      iter++;
    }
  }
  _hits.clear();
  _crystalhits.clear();
  _tracks.clear();
}

void DataInterface::removeAllComponents()
{
  _components.clear();
  _straws.clear();
  _crystals.clear();
  _hits.clear();
  _crystalhits.clear();
  _tracks.clear();
  _supportstructures.clear();
  _otherstructures.clear();
  if(_geometrymanager) delete _geometrymanager;
  _geometrymanager=NULL;
}

}
