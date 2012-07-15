#define USETRAJECTORY
#include "DataInterface.h"

#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Cube.h"
#include "Cylinder.h"
#include "EventDisplayFrame.h"
#include "FilterDialog.h"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "HepPID/ParticleName.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/PointTrajectoryCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "Straw.h"
#include "TTrackerGeom/inc/TTracker.hh"
#include "ITrackerGeom/inc/ITracker.hh"
#include "ITrackerGeom/inc/SuperLayer.hh"
#include "ITrackerGeom/inc/Cell.hh"
#include "ITrackerGeom/inc/CellGeometryHandle.hh"
#include "TargetGeom/inc/Target.hh"
#include "Track.h"
#include "TrackColorSelector.h"
#include "TrackerGeom/inc/Tracker.hh"
#include "art/Framework/Principal/Run.h"
#include "cetlib/map_vector.h"
#include "dict_classes/ComponentInfo.h"
#include "dict_classes/EventDisplayViewSetup.h"
#include <TAxis3D.h>
#include <TGFrame.h>
#include <TGeoVolume.h>
#include <TMath.h>
#include <TView.h>
#include <TGraphErrors.h>

#include <boost/shared_array.hpp>

#ifdef BABARINSTALLED
using namespace CLHEP;
#include "TrkBase/TrkRecoTrk.hh"
#include "TrkBase/TrkHotList.hh"
#include "KalmanTests/inc/TrkRecoTrkCollection.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "KalmanTrack/KalRep.hh"
#else
#warning BaBar package is absent. TrkRecoTrk cannot be displayed in the event display.
#endif

namespace mu2e_eventdisplay
{

DataInterface::DataInterface(EventDisplayFrame *mainframe):
              _geometrymanager(NULL),_topvolume(NULL),_mainframe(mainframe),
              _showUnhitStraws(false), _showUnhitCrystals(false)
{
    _minPoints=0;
    _minTime=NAN;
    _maxTime=NAN;
    _minMomentum=0;
    _showElectrons=true;
    _showMuons=true;
    _showGammas=true;
    _showNeutrinos=true;
    _showNeutrons=true;
    _showOthers=true;
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
    std::map<int,boost::shared_ptr<Straw> >::const_iterator straw;  //iterate over all straws
    for(straw=_straws.begin(); straw!=_straws.end(); straw++)
    {
      straw->second->update(time);
    }
  }
  else
  {
    std::vector<boost::shared_ptr<Straw> >::const_iterator hit; //iterate only over hit straws
    for(hit=_hits.begin(); hit!=_hits.end(); hit++)
    {
      (*hit)->update(time);
    }
  }
  if(_showUnhitCrystals)
  {
    std::map<int,boost::shared_ptr<Cube> >::const_iterator crystal; //iterate over all crystals
    for(crystal=_crystals.begin(); crystal!=_crystals.end(); crystal++) 
    {
      crystal->second->update(time);
    }
  }
  else
  {
    std::vector<boost::shared_ptr<Cube> >::const_iterator crystalhit; //iterate only over hit crystals
    for(crystalhit=_crystalhits.begin(); crystalhit!=_crystalhits.end(); crystalhit++)
    {
      (*crystalhit)->update(time);
    }
  }
  std::vector<boost::shared_ptr<Cylinder> >::const_iterator driftradius;
  for(driftradius=_driftradii.begin(); driftradius!=_driftradii.end(); driftradius++)
  {
    (*driftradius)->update(time);
  }
}

void DataInterface::setTrackFilter(unsigned int minPoints, double minTime, double maxTime, double minMomentum,
                                   bool showElectrons, bool showMuons, bool showGammas, 
                                   bool showNeutrinos, bool showNeutrons, bool showOthers)
{
    _minPoints=minPoints;
    _minTime=minTime;
    _maxTime=maxTime;
    _minMomentum=minMomentum;
    _showElectrons=showElectrons;
    _showMuons=showMuons;
    _showGammas=showGammas;
    _showNeutrinos=showNeutrinos;
    _showNeutrons=showNeutrons;
    _showOthers=showOthers;
}

void DataInterface::filterTracks()
{
  new FilterDialog(gClient->GetRoot(), this,
                   _minPoints, _minTime, _maxTime, _minMomentum,
                   _showElectrons, _showMuons, _showGammas, 
                   _showNeutrinos, _showNeutrons, _showOthers);
  std::vector<boost::shared_ptr<Track> >::const_iterator track;
  for(track=_tracks.begin(); track!=_tracks.end(); track++)
  {
    (*track)->setFilter(_minPoints, _minTime, _maxTime, _minMomentum,
                        _showElectrons, _showMuons, _showGammas, _showNeutrinos, _showNeutrons, _showOthers);
  }
}

void DataInterface::createGeometryManager()
{
  _geometrymanager = new TGeoManager("GeoManager", "GeoManager");
  _geometrymanager->SetVerboseLevel(0); 
  _topvolume = _geometrymanager->MakeBox("TopVolume", NULL, 1000, 1000, 1500);
  _geometrymanager->SetTopVolume(_topvolume);
  _geometrymanager->SetTopVisible(false);
  _geometrymanager->CloseGeometry();
  _geometrymanager->SetVisLevel(4);
  _topvolume->SetVisibility(0);
  _topvolume->SetLineColor(0);
  _topvolume->Draw("ogle");
  EventDisplayViewSetup::setup();
}

void DataInterface::fillGeometry()
{
  removeAllComponents();
  createGeometryManager();
  resetBoundaryP(_trackerMinmax);
  resetBoundaryP(_targetMinmax);
  resetBoundaryP(_calorimeterMinmax);
  resetBoundaryP(_tracksMinmax);

  art::ServiceHandle<mu2e::GeometryService> geom;

  const mu2e::SimpleConfig &config = geom->config();
  _xOffset=config.getDouble("mu2e.solenoidOffset");    //between Mu2e and Tracker coordinates
  _zOffset=-config.getDouble("mu2e.detectorSystemZ0"); //between Mu2e and Tracker coordinates
  _zOffsetDS=1800.0;                                   //between DS and Tracker coordinates

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
      int idStraw =  s.id().getStraw();
      int idLayer =  s.id().getLayer();
      int idSector =  s.id().getSector();
      int idDevice =  s.id().getDevice();
      int index = s.index().asInt();

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
    double innerRadius=ttracker->getSupportParams().innerRadius();
    double outerRadius=ttracker->getSupportParams().outerRadius();
    double zHalfLength=ttracker->getTrackerEnvelopeParams().zHalfLength();
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
                                          zHalfLength,innerRadius,outerRadius, NAN, 
                                          _geometrymanager, _topvolume, _mainframe, info, true));
    _components.push_back(shape);
    _supportstructures.push_back(shape);

//Envelope
    innerRadius=ttracker->getTrackerEnvelopeParams().innerRadius();
    outerRadius=ttracker->getTrackerEnvelopeParams().outerRadius();
    zHalfLength=ttracker->getTrackerEnvelopeParams().zHalfLength();

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
                                                  zHalfLength,innerRadius,outerRadius, NAN,
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
                                               zHalfLength,innerRadius,outerRadius, NAN, 
                                               _geometrymanager, _topvolume, _mainframe, infoToyDS, true));
    _components.push_back(shapeToyDS);
    _otherstructures.push_back(shapeToyDS);
  } else if(geom->hasElement<mu2e::ITracker>()) {
//Cells
    //mu2e::GeomHandle<mu2e::ITracker> itracker;

    const mu2e::Tracker& tracker = mu2e::getTrackerOrThrow();
    const mu2e::ITracker &itracker = static_cast<const mu2e::ITracker&>( tracker );

    double zoff = itracker.z0()-config.getDouble("mu2e.detectorSystemZ0");

    mu2e::CellGeometryHandle *itwp = itracker.getCellGeometryHandle();

    const boost::shared_array<mu2e::SuperLayer> sprlr=itracker.getSuperLayersArray();
    for ( int iS=0; iS<itracker.nSuperLayers(); iS++) {
            for ( int iL=0; iL<sprlr[iS].nLayers(); iL++ ) {
                    boost::shared_ptr<mu2e::ITLayer> ilr = sprlr[iS].getLayer(iL);
                    for ( int iC=0; iC<ilr->nCells(); iC++ ) {
                            mu2e::Cell *s = ilr->getCell(iC).get();

                            int idLayer =  s->Id().getLayer();
                            itwp->SelectCell(iS,idLayer,iC);

//                            const CLHEP::Hep3Vector& p = itwp->GetWireCenter(); //s->getMidPoint();
//                            const CLHEP::Hep3Vector& d = itwp->GetWireDirection(); //s->getDirection();
//                            double theta = d.theta();
//                            double phi = d.phi();

                            if (itracker.isDumbbell()){
                                    /*
                                    float wCntPos[3];
                                    double l = (itracker.zHalfLength()-itracker.zZonesLimits()[1])*0.5;
                                    float zWCnt = itracker.zZonesLimits()[1]+l;
                                    itwp->WirePosAtZ(zWCnt,wCntPos);
                                    wCntPos[2]+=zoff;
                                    double wl=l/cos(itwp->GetWireEpsilon()) ;
                                    */
                                    const CLHEP::Hep3Vector& p = itwp->GetCellCenter(); //s->getMidPoint();
                                    const CLHEP::Hep3Vector& d = itwp->GetCellDirection(); //s->getDirection();
                                    double theta = d.theta();
                                    double phi = d.phi();
                                    double x = p.x();
                                    double y = p.y();
                                    double z = p.z()+zoff;
                                    double l = itwp->GetCellHalfLength();
                                    int index = itwp->computeDet(iS,idLayer,iC); //s->index().asInt();

                                    char c[200];
                                    sprintf(c,"Cell %i  Layer %i  SuperLayer %i DownStream",iC,idLayer,iS);
                                    boost::shared_ptr<ComponentInfo> infoDnS(new ComponentInfo());
                                    infoDnS->setName(c);
                                    infoDnS->setText(0,c);
                                    //boost::shared_ptr<Straw> shapeDnS(new Straw(wCntPos[0],wCntPos[1],wCntPos[2], NAN, theta, phi, wl, _geometrymanager, _topvolume, _mainframe, infoDnS, false));
                                    boost::shared_ptr<Straw> shapeDnS(new Straw(x,y,z, NAN, theta, phi, l, _geometrymanager, _topvolume, _mainframe, infoDnS, false));
                                    _components.push_back(shapeDnS);
                                    _straws[index]=shapeDnS;

                                    /*
                                    zWCnt = itracker.zZonesLimits()[0]-l;
                                    itwp->WirePosAtZ(zWCnt,wCntPos);
                                    wCntPos[2]+=zoff;
                                    */
                                    itwp->SelectCell(iS,idLayer,iC,true);
                                    const CLHEP::Hep3Vector& pUp = itwp->GetCellCenter(); //s->getMidPoint();
                                    const CLHEP::Hep3Vector& dUp = itwp->GetCellDirection(); //s->getDirection();
                                    theta = dUp.theta();
                                    phi = dUp.phi();
                                    x = pUp.x();
                                    y = pUp.y();
                                    z = pUp.z()+zoff;
                                    l = itwp->GetCellHalfLength();
                                    index = itwp->computeDet(iS,idLayer,iC,true); //s->index().asInt();

                                    sprintf(c,"Cell %i  Layer %i  SuperLayer %i UpStream",iC,idLayer,iS);
                                    boost::shared_ptr<ComponentInfo> infoUpS(new ComponentInfo());
                                    infoUpS->setName(c);
                                    infoUpS->setText(0,c);
                                    //boost::shared_ptr<Straw> shapeUpS(new Straw(wCntPos[0],wCntPos[1],wCntPos[2], NAN, theta, phi, wl, _geometrymanager, _topvolume, _mainframe, infoUpS, false));
                                    boost::shared_ptr<Straw> shapeUpS(new Straw(x,y,z, NAN, theta, phi, l, _geometrymanager, _topvolume, _mainframe, infoUpS, false));
                                    _components.push_back(shapeUpS);
                                    _straws[index]=shapeUpS;

                            }else {
                                    const CLHEP::Hep3Vector& p = itwp->GetCellCenter(); //s->getMidPoint();
                                    const CLHEP::Hep3Vector& d = itwp->GetCellDirection(); //s->getDirection();
                                    double theta = d.theta();
                                    double phi = d.phi();
                                    double x = p.x();
                                    double y = p.y();
                                    double z = p.z()+zoff;
                                    double l = s->getHalfLength();
                                    int index = itwp->computeDet(iS,idLayer,iC); //s->index().asInt();

                                    char c[200];
                                    sprintf(c,"Cell %i  Layer %i  SuperLayer %i",iC,idLayer,iS);
                                    boost::shared_ptr<ComponentInfo> info(new ComponentInfo());
                                    info->setName(c);
                                    info->setText(0,c);
                                    boost::shared_ptr<Straw> shape(new Straw(x,y,z, NAN, theta, phi, l, _geometrymanager, _topvolume, _mainframe, info, false));
                                    _components.push_back(shape);
                                    _straws[index]=shape;
                            }
                    }
            }
    }
/*
//Support Structure
    double innerRadius=itracker.getSupportParams().innerRadius();
    double outerRadius=itracker.getSupportParams().outerRadius();
    double zHalfLength=itracker.getTrackerEnvelopeParams().zHalfLength();
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
                                          zHalfLength,innerRadius,outerRadius, NAN,
                                          _geometrymanager, _topvolume, _mainframe, info, true));
    _components.push_back(shape);
    _supportstructures.push_back(shape);
*/
//Envelope
    double innerRadius=itracker.r0();//getTrackerEnvelopeParams().innerRadius();
    double outerRadius=itracker.rOut();//getTrackerEnvelopeParams().outerRadius();
    double zHalfLength=itracker.maxEndCapDim();

    char c[200];
    boost::shared_ptr<ComponentInfo> infoEnvelope(new ComponentInfo());
    sprintf(c,"ITracker Envelope");
    infoEnvelope->setName(c);
    infoEnvelope->setText(0,c);
    sprintf(c,"Inner Radius %.f mm  Outer Radius %.f mm",innerRadius/CLHEP::mm,outerRadius/CLHEP::mm);
    infoEnvelope->setText(1,c);
    sprintf(c,"Length %.f mm",2.0*zHalfLength/CLHEP::mm);
    infoEnvelope->setText(2,c);
    sprintf(c,"Center at x: 0 mm, y: 0 mm, z: %f mm",zoff);
    infoEnvelope->setText(3,c);
    boost::shared_ptr<Cylinder> shapeEnvelope(new Cylinder(0,0,zoff, 0,0,0,
                                                  zHalfLength,innerRadius,outerRadius, NAN,
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
                                               zHalfLength,innerRadius,outerRadius, NAN,
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
    boost::shared_ptr<Cylinder> shape(new Cylinder(0,0,z, 0,0,0, length/2.0,0,radius, NAN, 
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
      int    id=v.id();
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
      sprintf(c,"Dimension  ?x: %.f mm, ?y: %.f mm, ?z: %.f mm",sx/CLHEP::mm,sy/CLHEP::mm,sz/CLHEP::mm);
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
      sprintf(c,"Dimension  ?x: %.f mm, ?y: %.f mm, ?z: %.f mm",sx/CLHEP::mm,crystalHalfSize/CLHEP::mm,crystalHalfSize/CLHEP::mm);
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
  //MBS
  if(config.getBool("hasMBS", false)) {
    char c[200];
    double mbsinr[3], mbsoutr[3], mbslen[3], mbsz[3];
    mbsinr[0]  = config.getDouble("mbs.BSTCInnerRadius");
    mbsoutr[0] = config.getDouble("mbs.BSTSOuterRadius");
    mbsinr[1]  = config.getDouble("mbs.BSBSInnerRadius");
    mbsoutr[1] = config.getDouble("mbs.BSTSOuterRadius");
    mbsinr[2]  = config.getDouble("mbs.CLV2InnerRadius");
    mbsoutr[2] = config.getDouble("mbs.CLV2OuterRadius");
    mbslen[0]  = config.getDouble("mbs.BSTCHLength");
    mbslen[1]  = config.getDouble("mbs.BSTSHLength") - mbslen[0];
    mbslen[2]  = config.getDouble("mbs.CLV2HLength");
    double mbsbstsz = config.getDouble("mbs.BSTSZ");
    double mbstotallen = (mbslen[0] + mbslen[1])*2.;
    double mbsstartz = mbsbstsz - mbstotallen/2. + _zOffset;
    mbsz[0] = mbsstartz + mbslen[0] ;
    mbsz[1] = mbsstartz + mbslen[0]*2. + mbslen[1];
    mbsz[2] = mbsstartz +mbstotallen - mbslen[2];
    for (unsigned int i = 0 ; i <3 ; ++i) {
    //for (unsigned int i = 2 ; i <3 ; ++i) {
      boost::shared_ptr<ComponentInfo> infoMBS(new ComponentInfo());
      sprintf(c,"MBS %d", i);
      infoMBS->setName(c);
      infoMBS->setText(0,c);
      sprintf(c,"Inner Radius %.f mm  Outer Radius %.f mm",mbsinr[i]/CLHEP::mm,mbsoutr[i]/CLHEP::mm);
      infoMBS->setText(1,c);
      sprintf(c,"Length %.f mm",2.0*mbslen[i]/CLHEP::mm);
      infoMBS->setText(2,c);
      sprintf(c,"Center at x: 0 mm, y: 0 mm, z: %.f mm",mbsz[i]/CLHEP::mm);
      infoMBS->setText(3,c);
      boost::shared_ptr<Cylinder> shapeMBS(new Cylinder(0,0,mbsz[i], 0,0,0,
                                               mbslen[i], mbsinr[i], mbsoutr[i], NAN,
                                               _geometrymanager, _topvolume, _mainframe, infoMBS, false));
      _components.push_back(shapeMBS);
      _mbsstructures.push_back(shapeMBS);
    }

  }

//   if(geom->hasElement<mu2e::CosmicRayShield>())
//   {
//     mu2e::GeomHandle<mu2e::CosmicRayShield> crs;

//     std::map<std::string,mu2e::CRSSteelShield> const & shields =
//       crs->getCRSSteelShields();

//     for (std::map<std::string,mu2e::CRSSteelShield>::const_iterator ishield=shields.begin();
//          ishield!=shields.end(); ++ishield) {
//       mu2e::CRSSteelShield const & steelshield = ishield->second;
//       std::string shieldName = ishield->first;
//       double x=steelshield.getGlobalOffset().x()-_mu2eOriginInWorld.x()+_xOffset;
//       double y=steelshield.getGlobalOffset().y()-_mu2eOriginInWorld.y();
//       double z=steelshield.getGlobalOffset().z()-_mu2eOriginInWorld.z()+_zOffset;
//       double dx=steelshield.getHalfLengths()[0];
//       double dy=steelshield.getHalfLengths()[1];
//       double dz=steelshield.getHalfLengths()[2];
//       char c[200];
//       boost::shared_ptr<ComponentInfo> info(new ComponentInfo());
//       sprintf(c,"Cosmic Ray Steel Shield (%s)",shieldName.c_str());
//       info->setName(c);
//       info->setText(0,c);
//       sprintf(c,"Dimension  ?x: %.f mm, ?y: %.f mm, ?z: %.f mm",dx/CLHEP::mm,dy/CLHEP::mm,dz/CLHEP::mm);
//       info->setText(1,c);
//       sprintf(c,"Center at x: %.f mm, y: %.f mm, z: %.f mm",x/CLHEP::mm,y/CLHEP::mm,z/CLHEP::mm);
//       info->setText(2,c);

//       double holeRadius=steelshield.getHoleRadius();
//       if(holeRadius==0)
//       {
//         boost::shared_ptr<Cube> shape(new Cube(x,y,z,  dx,dy,dz,  0,0,0, NAN,0,
//                                           _geometrymanager, _topvolume, _mainframe, info, true));
//         _components.push_back(shape);
//         _otherstructures.push_back(shape);
//       }
//       else
//       {
// //TODO: This needs to be replaced by a shape with a hole
//         boost::shared_ptr<Cube> shape(new Cube(x,y,z,  dx,dy,dz,  0,0,0, NAN,0,
//                                           _geometrymanager, _topvolume, _mainframe, info, true));
//         _components.push_back(shape);
//         _otherstructures.push_back(shape);
//       }
//     }
//   }

}

void DataInterface::makeMuonBeamStopStructuresVisible(bool visible)
{
  std::vector<boost::shared_ptr<VirtualShape> >::const_iterator structure;
  for(structure=_mbsstructures.begin(); structure!=_mbsstructures.end(); structure++)
  {
    (*structure)->setDefaultVisibility(visible);
    (*structure)->start();
  }

  //tracks and straws don't have to be pushed into the foreground if the structure is removed
  if(visible) toForeground();
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
  std::vector<boost::shared_ptr<Cylinder> >::const_iterator driftradius;
  for(driftradius=_driftradii.begin(); driftradius!=_driftradii.end(); driftradius++)
  {
    double time=(*driftradius)->getStartTime();
    if(hitcolors)
    {
      int color=TMath::FloorNint(20.0*(time-_hitsTimeMinmax.mint)/(_hitsTimeMinmax.maxt-_hitsTimeMinmax.mint));
      if(color>=20) color=19;
      if(color<=0) color=0;
      color+=2000;
      (*driftradius)->setColor(color);
    }
    else (*driftradius)->setColor(whitebackground?1:0);
  }
}

void DataInterface::useTrackColors(boost::shared_ptr<ContentSelector> const &contentSelector, bool trackcolors, bool whitebackground)
{
  std::vector<ContentSelector::trackInfoStruct> selectedTracks=contentSelector->getSelectedTrackNames();
  TrackColorSelector colorSelector(&selectedTracks);
  std::vector<boost::shared_ptr<Track> >::const_iterator track;
  for(track=_tracks.begin(); track!=_tracks.end(); track++)
  {
    if(trackcolors)
    {
      int color=colorSelector.getColor(*track);
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

void DataInterface::fillEvent(boost::shared_ptr<ContentSelector> const &contentSelector)
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


#ifdef BABARINSTALLED
  const mu2e::TrkRecoTrkCollection *trkRecoTrkHits=contentSelector->getSelectedHitCollection<mu2e::TrkRecoTrkCollection>();
  if(trkRecoTrkHits!=NULL)
  {
    boost::shared_ptr<TGraphErrors> residualGraph(new TGraphErrors());
    residualGraph->SetTitle("Residual Graph");

    for(unsigned int i=0; i<trkRecoTrkHits->size(); i++)
    {
      const TrkRecoTrk &particle = *trkRecoTrkHits->at(i);
      const TrkHotList* hots = particle.hots();
      if(hots!=NULL)
      {
        _numberHits+=hots->nHit();
        for(TrkHotList::hot_iterator iter=hots->begin(); iter!=hots->end(); iter++)
        {
          const TrkHitOnTrk *hitOnTrack = iter.get();
          const mu2e::TrkStrawHit* strawHit = dynamic_cast<const mu2e::TrkStrawHit*>(hitOnTrack);
          if(strawHit)
          {
            int    strawindex=strawHit->straw().index().asInt();
            double time = strawHit->time(); 
            double hitT0 = strawHit->hitT0(); //this is the time the hit "arrived at the straw"
                                              //don't know what the other times are
            double strawtime = strawHit->strawHit().time();
            double driftRadius = strawHit->driftRadius();
            const HepPoint &p=strawHit->hitTraj()->position(strawHit->hitLen());
            double theta = strawHit->straw().getDirection().theta();
            double phi = strawHit->straw().getDirection().phi();

            double residual, residualError;
            if(strawHit->resid(residual, residualError))
            {
              int n=residualGraph->GetN();
              residualGraph->SetPoint(n,p.z(),residual);
              residualGraph->SetPointError(n,0,residualError);
            }

            std::map<int,boost::shared_ptr<Straw> >::iterator straw=_straws.find(strawindex);
            if(straw!=_straws.end() && !isnan(time))
            {
              double previousStartTime=straw->second->getStartTime();
              if(isnan(previousStartTime))
              {
                findBoundaryT(_hitsTimeMinmax, hitT0);  //is it Ok to exclude all following hits from the time window?
                straw->second->setStartTime(hitT0);
                straw->second->start();
                char c1[100], c2[100], c3[100];
                sprintf(c1,"hitT0(s): %gns",hitT0/CLHEP::ns);
                sprintf(c2,"hit time(s): %gns",time/CLHEP::ns);
                sprintf(c3,"strawhit time(s): %gns",strawtime/CLHEP::ns);
                straw->second->getComponentInfo()->setText(1,c1);
                straw->second->getComponentInfo()->setText(2,c2);
                straw->second->getComponentInfo()->setText(3,c3);
                residualGraph->GetXaxis()->SetTitle("z [mm]");
                residualGraph->GetYaxis()->SetTitle("Residual [??]");
                straw->second->getComponentInfo()->getHistVector().push_back(boost::dynamic_pointer_cast<TObject>(residualGraph));
                _hits.push_back(straw->second);
              }
              else
              {
                straw->second->getComponentInfo()->expandLine(1,hitT0/CLHEP::ns,"ns");
                straw->second->getComponentInfo()->expandLine(2,time/CLHEP::ns,"ns");
                straw->second->getComponentInfo()->expandLine(3,strawtime/CLHEP::ns,"ns");
              }

              char c0[200], c1[200];
              const boost::shared_ptr<std::string> strawname=straw->second->getComponentInfo()->getName();
              sprintf(c0,"Drift Radius for %s",strawname->c_str());
              sprintf(c1,"Drift Radius %gcm",driftRadius/CLHEP::cm);
              boost::shared_ptr<ComponentInfo> info(new ComponentInfo());
              info->setName(c0);
              info->setText(0,strawname->c_str());
              info->setText(1,c1);
              boost::shared_ptr<Cylinder> driftradius(new Cylinder(p.x(),p.y(),p.z(), 
                                                          phi+TMath::Pi()/2.0,theta,0,
                                                          5, //the halflength of 5 has no meaning 
                                                          0,driftRadius,hitT0, 
                                                          _geometrymanager, _topvolume, _mainframe, info, false));
              _components.push_back(driftradius);
              _driftradii.push_back(driftradius);
            }
          }
        }
      }
    }
  }
#endif


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
  art::ServiceHandle<mu2e::GeometryService> geoservice;
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
  std::vector<ContentSelector::trackInfoStruct> trackInfos;
  std::vector<const mu2e::SimParticleCollection*> simParticleCollectionVector=contentSelector->getSelectedTrackCollection<mu2e::SimParticleCollection>(trackInfos);
  for(unsigned int i=0; i<simParticleCollectionVector.size(); i++)
  {
    const mu2e::SimParticleCollection *simParticles=simParticleCollectionVector[i];
    cet::map_vector<mu2e::SimParticle>::const_iterator iter;
    for(iter=simParticles->begin(); iter!=simParticles->end(); iter++)
    {
      const mu2e::SimParticle& particle = iter->second;
      const cet::map_vector_key& particleKey = iter->first;
      int id = particle.id().asInt();
      int parentid = -1;
      if(particle.hasParent()) parentid = particle.parent()->id().asInt();
      int particleid=particle.pdgId();
      int trackclass=trackInfos[i].classID;
      int trackclassindex=trackInfos[i].index;
      std::string particlecollection=trackInfos[i].entryText;
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
      sprintf(c0,"Track %i  %s  (%s)",id,particlename.c_str(),particlecollection.c_str());
      if(parentid!=0) sprintf(c1,"Track %i  %s  Parent %i  (%s)",id,particlename.c_str(),parentid,particlecollection.c_str());
      else sprintf(c1,"Track %i  %s  generated track  (%s)",id,particlename.c_str(),particlecollection.c_str());
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
      std::vector<cet::map_vector_key>::const_iterator daughter;
      for(daughter=particle.daughterIds().begin();
          daughter!=particle.daughterIds().end();
          daughter++)
      {
        info->expandLine(4,daughter->asInt(),"");
      }
      boost::shared_ptr<Track> shape(new Track(x1,y1,z1,t1, x2,y2,z2,t2, particleid, trackclass, trackclassindex, e1, 
                                               _geometrymanager, _topvolume, _mainframe, info));
      findTrajectory(contentSelector,shape,particleKey, t1,t2, simParticles,particle.daughterIds(), trackInfos[i]);
      _components.push_back(shape);
      _tracks.push_back(shape);
    }
  }

#ifdef BABARINSTALLED
  trackInfos.clear();
  std::vector<const mu2e::TrkRecoTrkCollection*> trkRecoTrkCollectionVector=contentSelector->getSelectedTrackCollection<mu2e::TrkRecoTrkCollection>(trackInfos);
  for(unsigned int i=0; i<trkRecoTrkCollectionVector.size(); i++)
  {
    const mu2e::TrkRecoTrkCollection *trkRecoTrks=trkRecoTrkCollectionVector[i];
    for(unsigned int j=0; j<trkRecoTrks->size(); j++)
    {
      const TrkRecoTrk &particle = *trkRecoTrks->at(j);
      double t0=particle.trackT0();
      const TrkRep* trkrep = particle.getRep(particle.defaultType());
      if(trkrep!=NULL)
      {
        int particleid=0;
        int trackclass=trackInfos[i].classID;
        int trackclassindex=trackInfos[i].index;
        std::string trackcollection=trackInfos[i].entryText;
        switch(static_cast<int>(trkrep->particleType()))
        {
           case 0 : particleid=11;   break;  //electron
           case 1 : particleid=13;   break;  //muon
           case 2 : particleid=211;  break;  //pion
           case 3 : particleid=321;  break;  //kaon
           case 4 : particleid=2212; break;  //proton
        };
        std::string particlename=HepPID::particleName(particleid);
        char c0[200], c2[200], c3[200], c4[200];
        sprintf(c0,"Kalman Track %i  %s  (%s)",static_cast<int>(particle.id()),particlename.c_str(),trackcollection.c_str());
        boost::shared_ptr<ComponentInfo> info(new ComponentInfo());
        info->setName(c0);
        info->setText(0,c0);

        double hitcount=0;
        double offset=0;
        const TrkHotList* hots=trkrep->hotList();
        if(hots!=NULL)
        {
          boost::shared_ptr<TGraphErrors> residualGraph(new TGraphErrors());
          residualGraph->SetTitle("Residual Graph");
          for(TrkHotList::hot_iterator iter=hots->begin(); iter!=hots->end(); iter++)
          {
            const TrkHitOnTrk *hitOnTrack = iter.get();
            const mu2e::TrkStrawHit* strawHit = dynamic_cast<const mu2e::TrkStrawHit*>(hitOnTrack);
            if(strawHit)
            {
              double strawTime   = strawHit->hitT0()/CLHEP::ns;
              double trackTime   = strawTime*CLHEP::ns;  //TODO: add correction for drift time
              double weight= strawHit->weight();  
              double fltLen= strawHit->fltLen();
              const HepPoint &p=strawHit->hitTraj()->position(strawHit->hitLen());
              double t     = trkrep->arrivalTime(fltLen);
              offset+=(trackTime-t)*weight;
              hitcount+=weight;

              double residual, residualError;
              if(strawHit->resid(residual, residualError))
              {
                int n=residualGraph->GetN();
                residualGraph->SetPoint(n,p.z(),residual);
                residualGraph->SetPointError(n,0,residualError);
              }
            }
          }
          residualGraph->GetXaxis()->SetTitle("z [mm]");
          residualGraph->GetYaxis()->SetTitle("Residual [??]");
          info->getHistVector().push_back(boost::dynamic_pointer_cast<TObject>(residualGraph));
        }
        if(hitcount>0) offset/=hitcount; else offset=0;

        double fltLMin=trkrep->startValidRange();
        double fltLMax=trkrep->endValidRange();
        double p1=trkrep->momentum(fltLMin).mag();
        double p2=trkrep->momentum(fltLMax).mag();
        double x1=trkrep->position(fltLMin).x();
        double y1=trkrep->position(fltLMin).y();
        double z1=trkrep->position(fltLMin).z();
        double x2=trkrep->position(fltLMax).x();
        double y2=trkrep->position(fltLMax).y();
        double z2=trkrep->position(fltLMax).z();
        double t1=trkrep->arrivalTime(fltLMin)+offset;
        double t2=trkrep->arrivalTime(fltLMax)+offset;
        boost::shared_ptr<Track> track(new Track(x1,y1,z1,t1, x2,y2,z2,t2, particleid, trackclass, trackclassindex, p1, 
                                                 _geometrymanager, _topvolume, _mainframe, info));
        _components.push_back(track);
        _tracks.push_back(track);

        double fltStep = (fltLMax - fltLMin)/400.0;
        for(unsigned int step = 0; step <= 400.0; step++) 
        {		
          double fltL = fltLMin + step*fltStep;
          double   t = trkrep->arrivalTime(fltL)+offset;
          HepPoint p = trkrep->position(fltL);
          findBoundaryT(_tracksTimeMinmax, t);
          findBoundaryP(_tracksMinmax, p.x(), p.y(), p.z());
          track->addTrajectoryPoint(p.x(), p.y(), p.z(), t);
        }
        const KalRep* kalrep = dynamic_cast<const KalRep*>(trkrep);
        if(kalrep!=NULL)
        {
          int charge = kalrep->charge();
          sprintf(c2,"Charge %i",charge);
          info->setText(1,c2);
          sprintf(c3,"Start Momentum %gMeV/c  End Momentum %gMeV/c",p1/CLHEP::MeV,p2/CLHEP::MeV);
          sprintf(c4,"T0 %gns",t0/CLHEP::ns);
          info->setText(2,c3);
          info->setText(3,c4);
        }
      }
    }
  }
#endif
}

void DataInterface::findTrajectory(boost::shared_ptr<ContentSelector> const &contentSelector,
                                   boost::shared_ptr<Track> const &track, const cet::map_vector_key &id,
                                   double t1, double t2,
                                   const mu2e::SimParticleCollection *simParticles,
                                   const std::vector<cet::map_vector_key> &daughterVect,
                                   const ContentSelector::trackInfoStruct &trackInfo)
{
#ifdef USETRAJECTORY
  const mu2e::PointTrajectoryCollection *pointTrajectories=contentSelector->getPointTrajectoryCollection(trackInfo);
  if(pointTrajectories!=NULL)
  {
    const mu2e::PointTrajectory* trajectory=pointTrajectories->getOrNull(id);
    if(trajectory)
    {
// fill a vector with all trajectory points, leave the times as NAN,
// except the first and last point, which get the times from the track
        std::vector<trajectoryStruct> trajectoryVect;
        const std::vector<CLHEP::Hep3Vector>& pVect=trajectory->points();
        for(unsigned int i=0; i<pVect.size(); i++)
        {
          trajectoryStruct ts;
          ts.v=pVect[i];
          ts.v.setZ(pVect[i].z()+_zOffsetDS);
          if(i==0) ts.t=t1;
          if(i==(pVect.size()-1)) ts.t=t2;
          trajectoryVect.push_back(ts);
        }

// the next section will disappear once the time information is present in the PointTrajectory class
// loop over all daughter ids
        std::vector<cet::map_vector_key>::const_iterator daughterIter;
        for(daughterIter=daughterVect.begin(); daughterIter!=daughterVect.end(); daughterIter++)
        {
// find the track which corresponds to each daughter id
          const mu2e::SimParticle *particle=simParticles->getOrNull(*daughterIter);
          if(particle)
          {
// extract the start position and start time of the daughter track
              CLHEP::Hep3Vector v;
              v.setX(particle->startPosition().x()+_xOffset);
              v.setY(particle->startPosition().y());
              v.setZ(particle->startPosition().z()+_zOffset);
              double newTime=particle->startGlobalTime();
              if(newTime>t2 || newTime<t1) break;
// try to find the point in the trajectory vector, which is closest to this starting point
              double mindiff=NAN;
              unsigned int mindiffPoint=0;
              for(unsigned int j=0; j<trajectoryVect.size(); j++)
              {
                CLHEP::Hep3Vector diff=trajectoryVect[j].v-v;
                double diff_magnitude=diff.getR();
                if(diff_magnitude<mindiff || isnan(mindiff))
                {
                  mindiff=diff_magnitude;
                  mindiffPoint=j;
                }
              }
// if a valid trajectory point is found, set the time of this trajectory point to the starting time of the daughter track
              if(mindiffPoint>0 && mindiffPoint<(trajectoryVect.size()-1))
              {
                double oldTime=trajectoryVect[mindiffPoint].t;
                if(newTime<oldTime || isnan(oldTime)) trajectoryVect[mindiffPoint].t=newTime;
              }
          }
        }

// loop over all trajectory points (some of them have their times set from the previous steps,
// the remaining points have their times equal to NAN)
        int lastTimeEntry=0;
        double lastTime=t1;
        for(unsigned int i=1; i<trajectoryVect.size(); i++)
        {
// find a point with a non-NAN time
          double nextTime=trajectoryVect[i].t;
          if(!isnan(nextTime))
          {
            if(nextTime<lastTime) {trajectoryVect[i].t=NAN; continue;}
            double timeIntervall=(nextTime-lastTime)/(i-lastTimeEntry);
            double timeStep=lastTime;
// fill all times between the previous point with a non-NAN time to the current point with a non-NAN time
            for(unsigned int j=lastTimeEntry+1; j<i; j++)
            {
              timeStep+=timeIntervall;
              trajectoryVect[j].t=timeStep;
            }
            lastTimeEntry=i;
            lastTime=nextTime;
          }
        }

//replace the track, which previously had only two points, with a complete trajectory
        for(unsigned int i=0; i<trajectoryVect.size(); i++)
        {
          track->addTrajectoryPoint(trajectoryVect[i].v.getX(),
                                    trajectoryVect[i].v.getY(),
                                    trajectoryVect[i].v.getZ(),
                                    trajectoryVect[i].t);
        }
    }
  }
#endif
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

  std::vector<boost::shared_ptr<Straw> >::const_iterator hit;
  for(hit=_hits.begin(); hit!=_hits.end(); hit++)
  {
    for(int i=1; i<5; i++) (*hit)->getComponentInfo()->removeLine(i);  //keep first line
    (*hit)->getComponentInfo()->getHistVector().clear();
  }
  std::vector<boost::shared_ptr<Cube> >::const_iterator crystalhit;
  for(crystalhit=_crystalhits.begin(); crystalhit!=_crystalhits.end(); crystalhit++)
  {
    for(int i=1; i<5; i++) (*crystalhit)->getComponentInfo()->removeLine(i); //keep first line
    (*crystalhit)->getComponentInfo()->getHistVector().clear();
  }

  _hits.clear();
  _crystalhits.clear();
  _tracks.clear();
  _driftradii.clear();

  _mainframe->getHistDrawVector().clear();
}

void DataInterface::removeAllComponents()
{
  _components.clear();
  _straws.clear();
  _crystals.clear();
  _hits.clear();
  _crystalhits.clear();
  _tracks.clear();
  _driftradii.clear();
  _supportstructures.clear();
  _otherstructures.clear();
  _mbsstructures.clear();
  delete _geometrymanager;
  _geometrymanager=NULL;

  _mainframe->getHistDrawVector().clear();
}

}
