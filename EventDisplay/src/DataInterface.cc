#include "DataInterface.h"
#include "Track.h"
#include "Straw.h"
#include "SupportTTracker.h"
#include "dict_classes/ComponentInfo.h"

#include <TView.h>
#include <TAxis3D.h>
#include <TGFrame.h>
#include <TMath.h>

#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "LTrackerGeom/inc/LTracker.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "ITrackerGeom/inc/ITracker.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "HepPID/ParticleName.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "ToyDP/inc/ToyGenParticleCollection.hh"
#include "ToyDP/inc/SimParticleCollection.hh"
#include "ToyDP/inc/PointTrajectoryCollection.hh"

#include <TGeoTrack.h>
#include <TGeoVolume.h>

namespace mu2e_eventdisplay
{

DataInterface::DataInterface(const TGMainFrame *mainframe):
              _geometrymanager(NULL),_topvolume(NULL),_mainframe(mainframe),
              _showUnhitStraws(false)
{
}

DataInterface::~DataInterface()
{
  removeAllComponents();
}

const std::list<boost::shared_ptr<VirtualShape> > &DataInterface::getComponents() 
{
  return _components;
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
  gPad->GetView()->SetParallel();
  int irep=0;
  gPad->GetView()->SetView(200,70,90,irep);
  gPad->SetPhi(-90-180);
  gPad->SetTheta(90-70);
  gPad->GetView()->ShowAxis();
  TAxis3D::GetPadAxis(gPad)->SetLabelSize(0.025); 
  gPad->GetView()->Draw();
  gPad->Modified();
  gPad->Update();
}

void DataInterface::fillGeometry()
{
  removeAllComponents();
  createGeometryManager();
  _hitsMinmax.minx=NAN;
  _hitsMinmax.miny=NAN;
  _hitsMinmax.minz=NAN;
  _hitsMinmax.maxx=NAN;
  _hitsMinmax.maxy=NAN;
  _hitsMinmax.maxz=NAN;

  edm::Service<mu2e::GeometryService> geom;
  if(geom->hasElement<mu2e::TTracker>())
  {
    mu2e::GeomHandle<mu2e::TTracker> ttracker;  
    _z0 = ttracker->z0();
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

      findBoundaryP(_hitsMinmax, x, y, z);

      char c[200];
      sprintf(c,"Straw %i  Layer %i  Sector %i  Device %i",idStraw,idLayer,idSector,idDevice);
      boost::shared_ptr<ComponentInfo> info(new ComponentInfo());
      info->setText(0,c);
      boost::shared_ptr<Straw> shape(new Straw(x,y,z, NAN, theta, phi, l, _geometrymanager, _topvolume, _mainframe, info, false));
      _components.push_back(shape);
      _straws[index]=shape;
    }

    char c[200];
    boost::shared_ptr<ComponentInfo> info(new ComponentInfo());
    sprintf(c,"TTracker Support Structure");
    info->setText(0,c);
    sprintf(c,"Inner Radius %g mm",ttracker->getSupportParams().innerRadius);
    info->setText(1,c);
    sprintf(c,"Outer Radius %g mm",ttracker->getSupportParams().outerRadius);
    info->setText(2,c);
    sprintf(c,"Length %g mm",2.0*ttracker->getTrackerEnvelopeParams().zHalfLength);
    info->setText(3,c);
    boost::shared_ptr<SupportTTracker> shape(new SupportTTracker(0,0,0, 0,0, 
                                   ttracker->getTrackerEnvelopeParams().zHalfLength,
                                   ttracker->getSupportParams().innerRadius,
                                   ttracker->getSupportParams().outerRadius,
                                   _geometrymanager, _topvolume, _mainframe, info, true));
    _components.push_back(shape);
    _supportstructures.push_back(shape);
  }

  if(geom->hasElement<mu2e::ITracker>())
  {
    mu2e::GeomHandle<mu2e::ITracker> itracker;  
    _z0 = itracker->z0();
//TODO
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

void DataInterface::useHitColors(bool hitcolors, bool whitebackground)
{
  std::vector<boost::shared_ptr<Straw> >::const_iterator hit;
  for(hit=_hits.begin(); hit!=_hits.end(); hit++)
  {
    double time=(*hit)->getStartTime();
    if(hitcolors)
    {
      int color=TMath::FloorNint(20.0*(time-_hitsMinmax.mint)/(_hitsMinmax.maxt-_hitsMinmax.mint));
      if(color>=20) color=19;
      if(color<=0) color=0;
      color+=2000;
      (*hit)->setColor(color);
    }
    else (*hit)->setColor(whitebackground?1:0);
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

void DataInterface::findBoundaryT(minmax &m, double t)
{
  if(isnan(m.mint) || t<m.mint) m.mint=t;
  if(isnan(m.maxt) || t>m.maxt) m.maxt=t;
}

void DataInterface::findBoundaryP(minmax &m, double x, double y, double z)
{
  if(isnan(m.minx) || x<m.minx) m.minx=x;
  if(isnan(m.miny) || y<m.miny) m.miny=y;
  if(isnan(m.minz) || z<m.minz) m.minz=z;
  if(isnan(m.maxx) || x>m.maxx) m.maxx=x;
  if(isnan(m.maxy) || y>m.maxy) m.maxy=y;
  if(isnan(m.maxz) || z>m.maxz) m.maxz=z;
}

void DataInterface::fillEvent(const edm::Event& event)
{
  removeNonGeometryComponents();
  if(!_geometrymanager) createGeometryManager();

  edm::Handle<mu2e::StepPointMCCollection> hits;
  std::string _g4ModuleLabel = "g4run";
  std::string _trackerStepPoints = "tracker"; //TODO: this may not always be correct 
                                         //in the future: let user decide via display
  if(event.getByLabel(_g4ModuleLabel,_trackerStepPoints,hits))  
                                         //TODO: is this return bool to be used like this?
  {
    _numberHits=hits->size();
    _hitsMinmax.mint=NAN;
    _hitsMinmax.maxt=NAN;
    std::vector<mu2e::StepPointMC>::const_iterator iter;
    for(iter=hits->begin(); iter!=hits->end(); iter++)
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
          findBoundaryT(_hitsMinmax, time);  //is it Ok to exclude all following hits from the time window?
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


  edm::Handle<mu2e::ToyGenParticleCollection> genParticles;
  if(event.getByType(genParticles));
  {
    for(unsigned int i=0; i<genParticles->size(); i++)
    {
      const mu2e::ToyGenParticle& genparticle = (*genParticles)[i];
      int particleid=genparticle.pdgId();
      std::string particlename=HepPID::particleName(genparticle.pdgId());
//TODO: figure out coordinate transformation
//can these numbers be extracted from somewhere?
//if not, put them somehwere else together with all other numbers
      double x1=genparticle._position.x()+3904.0;
      double y1=genparticle._position.y();
      double z1=genparticle._position.z()-10200.0;  //-z0; //-12000.0;
      double t1=genparticle._time;
      double e1=genparticle._momentum.e();
      double x2=NAN;
      double y2=NAN;
      double z2=NAN;
      double t2=NAN;
      double length=NAN;
      std::string daughterString="Daughter ID(s):";

      edm::Handle<mu2e::SimParticleCollection> simParticles;
      if(event.getByType(simParticles))
      {
        MapVector<mu2e::SimParticle>::const_iterator iter;
        for(iter=simParticles->begin(); iter!=simParticles->end(); iter++)
        {
          const mu2e::SimParticle& simparticle = iter->second;
          unsigned int genindex = simparticle.generatorIndex();
          if(genindex==i)
          {
            int id = simparticle.id().asInt();
            char c[10];
            if(isnan(length)) sprintf(c," %i",id);
            else sprintf(c,", %i",id);
            daughterString.append(c); 

            CLHEP::Hep3Vector d;
            d=simparticle.startPosition()-genparticle._position;
            if(length<d.mag() || isnan(length))
            {
//TODO: figure out coordinate transformation
//can these numbers be extracted from somewhere?
//if not, put them somehwere else together with all other numbers
              x2=simparticle.startPosition().x()+3904.0;
              y2=simparticle.startPosition().y();
              z2=simparticle.startPosition().z()-10200.0;  //-z0; //-12000.0;
              t2=simparticle.startGlobalTime();
              length=d.mag();
            }
          }
        }
      }

      if(!isnan(length))
      {
        char c1[200],c2[200],c3[200];
        sprintf(c1,"Generated Track %s",particlename.c_str());
        sprintf(c2,"Start Energy %gMeV",e1/CLHEP::MeV);
        sprintf(c3,"Track Length %gmm",length/CLHEP::mm);
        boost::shared_ptr<ComponentInfo> info(new ComponentInfo());
        info->setText(0,c1);
        info->setText(1,c2);
        info->setText(2,c3);
        info->setText(3,daughterString.c_str());
        boost::shared_ptr<Track> shape(new Track(x1,y1,z1,t1, x2,y2,z2,t2, particleid, _geometrymanager, _topvolume, _mainframe, info));
        _components.push_back(shape);
        _tracks.push_back(shape);
      }
    }
  }

  _tracksMinmax.minx=NAN;
  _tracksMinmax.miny=NAN;
  _tracksMinmax.minz=NAN;
  _tracksMinmax.mint=NAN;
  _tracksMinmax.maxx=NAN;
  _tracksMinmax.maxy=NAN;
  _tracksMinmax.maxz=NAN;
  _tracksMinmax.maxt=NAN;
  edm::Handle<mu2e::SimParticleCollection> simParticles;
  if(event.getByType(simParticles))
  {
    MapVector<mu2e::SimParticle>::const_iterator iter;
    for(iter=simParticles->begin(); iter!=simParticles->end(); iter++)
    {
      const mu2e::SimParticle& particle = iter->second;
      int id = particle.id().asInt();
      int parentid = particle.parentId().asInt();
      int particleid=particle.pdgId();
      std::string particlename=HepPID::particleName(particle.pdgId());
//TODO: figure out coordinate transformation
//can these numbers be extracted from somewhere?
//if not, put them somehwere else together with all other numbers
      double x1=particle.startPosition().x()+3904.0;
      double y1=particle.startPosition().y();
      double z1=particle.startPosition().z()-10200.0;  //-z0; //-12000.0;
      double t1=particle.startGlobalTime();
      double e1=particle.startMomentum().e();
      double x2=particle.endPosition().x()+3904.0;
      double y2=particle.endPosition().y();
      double z2=particle.endPosition().z()-10200.0; //-z0; //-12000.0;
      double t2=particle.endGlobalTime();
      double e2=particle.endMomentum().e();

      findBoundaryT(_tracksMinmax, t1);
      findBoundaryT(_tracksMinmax, t2);
      findBoundaryP(_tracksMinmax, x1, y1, z1);
      findBoundaryP(_tracksMinmax, x2, y2, z2);

      CLHEP::Hep3Vector d = particle.endPosition() - particle.startPosition();
      double length=d.mag();

      char c1[200],c2[200],c3[200];
      sprintf(c1,"Track %i  %s  Parent %i",id,particlename.c_str(),parentid);
      sprintf(c2,"Start Energy %gMeV  End Energy %gMeV",e1/CLHEP::MeV,e2/CLHEP::MeV);
      sprintf(c3,"Track Length %gmm",length/CLHEP::mm);
      boost::shared_ptr<ComponentInfo> info(new ComponentInfo());
      info->setText(0,c1);
      info->setText(1,c2);
      info->setText(2,c3);
      info->setText(3,"Daughter IDs:");
      std::vector<MapVectorKey>::const_iterator daughter;
      for(daughter=particle.daughterIds().begin(); 
          daughter!=particle.daughterIds().end();
          daughter++)
      {
        info->expandLine(3,daughter->asInt(),"");
      }
      boost::shared_ptr<Track> shape(new Track(x1,y1,z1,t1, x2,y2,z2,t2, particleid, _geometrymanager, _topvolume, _mainframe, info));
      findTrajectory(event,shape,id);
      _components.push_back(shape);
      _tracks.push_back(shape);
    }
  }
}

bool DataInterface::findTrajectory(const edm::Event& event, 
                                   boost::shared_ptr<Track> track, int id)
{
  edm::Handle<mu2e::PointTrajectoryCollection> pointTrajectories;
  if(event.getByType(pointTrajectories))
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
//TODO: figure out coordinate transformation
//can this number be extracted from somewhere?
//if not, put it somehwere else together with all other numbers
          track->addTrajectoryPoint(p_vec[i].x(),p_vec[i].y(),p_vec[i].z()+1800,i,p_vec.size());
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
  _tracks.clear();
}

void DataInterface::removeAllComponents()
{
  std::list<boost::shared_ptr<VirtualShape> >::iterator iter;
  _components.clear();
  _straws.clear();
  _hits.clear();
  _tracks.clear();
  _supportstructures.clear();
  if(_geometrymanager) delete _geometrymanager;
  _geometrymanager=NULL;
}

}
