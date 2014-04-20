#define USETRAJECTORY
//
// $Id: DataInterface.cc,v 1.73 2014/04/20 10:57:09 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2014/04/20 10:57:09 $
//

#include "DataInterface.h"

#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "CalorimeterGeom/inc/VaneCalorimeter.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "EventDisplay/src/Cube.h"
#include "EventDisplay/src/Cylinder.h"
#include "EventDisplay/src/Cone.h"
#include "EventDisplay/src/EventDisplayFrame.h"
#include "EventDisplay/src/Hexagon.h"
#include "EventDisplay/src/Straw.h"
#include "EventDisplay/src/Track.h"
#include "EventDisplay/src/TrackColorSelector.h"
#include "EventDisplay/src/dict_classes/ComponentInfo.h"
#include "EventDisplay/src/dict_classes/EventDisplayViewSetup.h"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "HepPID/ParticleName.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoMultiCollection.hh"
#include "MCDataProducts/inc/MCTrajectoryCollection.hh"
#include "MCDataProducts/inc/PointTrajectoryCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "Mu2eUtilities/inc/PhysicalVolumeMultiHelper.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/TrkExtTrajCollection.hh"
#include "RecoDataProducts/inc/TrkExtTraj.hh"
#include "RecoDataProducts/inc/TrkExtTrajCollection.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "ITrackerGeom/inc/ITracker.hh"
#include "ITrackerGeom/inc/SuperLayer.hh"
#include "ITrackerGeom/inc/Cell.hh"
#include "ITrackerGeom/inc/CellGeometryHandle.hh"
#include "StoppingTargetGeom/inc/StoppingTarget.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "art/Framework/Principal/Run.h"
#include "cetlib/map_vector.h"
#include <TAxis3D.h>
#include <TGFrame.h>
#include <TGeoVolume.h>
#include <TMath.h>
#include <TView.h>
#include <TGraphErrors.h>

#include <boost/shared_array.hpp>

#ifdef BABARINSTALLED
using namespace CLHEP;
#include "TrkBase/TrkHotList.hh"
#include "KalmanTests/inc/KalRepCollection.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "KalmanTrack/KalRep.hh"
#else
#warning BaBar package is absent. KalRep cannot be displayed in the event display.
#endif

namespace mu2e_eventdisplay
{

DataInterface::DataInterface(EventDisplayFrame *mainframe):
              _geometrymanager(nullptr),_topvolume(nullptr),_mainframe(mainframe)
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
//  removeNonGeometryComponents();
}

void DataInterface::updateComponents(double time, boost::shared_ptr<ContentSelector> contentSelector)
{
  std::vector<boost::shared_ptr<Track> >::const_iterator track;
  for(track=_tracks.begin(); track!=_tracks.end(); track++)
  {
    (*track)->setFilter(_minPoints, _minTime, _maxTime, _minMomentum,
                        _showElectrons, _showMuons, _showGammas, _showNeutrinos, _showNeutrons, _showOthers);
    (*track)->update(time);
  }

  const mu2e::StrawHitFlagCollection *hitFlagCollection=contentSelector->getStrawHitFlagCollection();
  std::vector<boost::shared_ptr<Straw> >::const_iterator hit; //iterate only over hit straws
  for(hit=_hits.begin(); hit!=_hits.end(); hit++)
  {
    int hitnumber=(*hit)->getHitNumber();
    if(hitFlagCollection!=nullptr)
    {
      //flag collection selected --> show only hits which have certain hit flags
      (*hit)->setFilter(_minTime, _maxTime, true);
      if(hitnumber<static_cast<int>(hitFlagCollection->size()) && hitnumber>=0)
      {
        const mu2e::StrawHitFlag& hitFlag = (*hitFlagCollection)[hitnumber];
        if(hitFlag.hasAnyProperty(_hitFlagSetting))
        {
          (*hit)->setFilter(_minTime, _maxTime, false);
        } 
      }
    }
    else
    {
      //no flag collection selected --> show all hits
      (*hit)->setFilter(_minTime, _maxTime, false);
    }
    (*hit)->update(time);
  }

  std::vector<boost::shared_ptr<VirtualShape> >::const_iterator crystalhit; //iterate only over hit crystals
  for(crystalhit=_crystalhits.begin(); crystalhit!=_crystalhits.end(); crystalhit++)
  {
    (*crystalhit)->setFilter(_minTime, _maxTime);
    (*crystalhit)->update(time);
  }
  
  std::vector<boost::shared_ptr<Cylinder> >::const_iterator driftradius;
  for(driftradius=_driftradii.begin(); driftradius!=_driftradii.end(); driftradius++)
  {
    (*driftradius)->setFilter(_minTime, _maxTime);
    (*driftradius)->update(time);
  }
}

void DataInterface::getFilterValues(unsigned int &minPoints, double &minTime, double &maxTime, double &minMomentum,
                                    bool &showElectrons, bool &showMuons, bool &showGammas, 
                                    bool &showNeutrinos, bool &showNeutrons, bool &showOthers,
                                    mu2e::StrawHitFlag &hitFlagSetting)
{
    minPoints=_minPoints;
    minTime=_minTime;
    maxTime=_maxTime;
    minMomentum=_minMomentum;
    showElectrons=_showElectrons;
    showMuons=_showMuons;
    showGammas=_showGammas;
    showNeutrinos=_showNeutrinos;
    showNeutrons=_showNeutrons;
    showOthers=_showOthers;
    hitFlagSetting=_hitFlagSetting;
}

void DataInterface::setFilterValues(unsigned int minPoints, double minTime, double maxTime, double minMomentum,
                                    bool showElectrons, bool showMuons, bool showGammas, 
                                    bool showNeutrinos, bool showNeutrons, bool showOthers,
                                    mu2e::StrawHitFlag hitFlagSetting)
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
    _hitFlagSetting=hitFlagSetting;
}

DataInterface::timeminmax DataInterface::getHitsTimeBoundary() 
{
  DataInterface::timeminmax toReturn=_hitsTimeMinmax;
  if(_minTime>toReturn.mint) toReturn.mint=_minTime;
  if(_maxTime<toReturn.maxt) toReturn.maxt=_maxTime;
  return toReturn;
}

DataInterface::timeminmax DataInterface::getTracksTimeBoundary() 
{
  DataInterface::timeminmax toReturn=_tracksTimeMinmax;
  if(_minTime>toReturn.mint) toReturn.mint=_minTime;
  if(_maxTime<toReturn.maxt) toReturn.maxt=_maxTime;
  return toReturn;
}

void DataInterface::createGeometryManager()
{
  _geometrymanager = new TGeoManager("GeoManager", "GeoManager");
  _geometrymanager->SetVerboseLevel(0); 
  _topvolume = _geometrymanager->MakeBox("TopVolume", nullptr, 1000, 1000, 1500);
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
  _detSysOrigin = mu2e::GeomHandle<mu2e::DetectorSystem>()->getOrigin();

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
      boost::shared_ptr<Straw> shape(new Straw(x,y,z, NAN, theta, phi, l, 
                                               _geometrymanager, _topvolume, _mainframe, info, true));
      _components.push_back(shape);
      _straws[index]=shape;
    }

//Support Structure
    double innerRadius=ttracker->getSupportParams().innerRadius();
    double outerRadius=ttracker->getSupportParams().outerRadius();
    double zHalfLength=ttracker->getInnerTrackerEnvelopeParams().zHalfLength();
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
    shape->makeGeometryVisible(true);
    _components.push_back(shape);
    _supportstructures.push_back(shape);

//Envelope
    innerRadius=ttracker->getInnerTrackerEnvelopeParams().innerRadius();
    outerRadius=ttracker->getInnerTrackerEnvelopeParams().outerRadius();
    zHalfLength=ttracker->getInnerTrackerEnvelopeParams().zHalfLength();

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
    shapeEnvelope->makeGeometryVisible(true);
    _components.push_back(shapeEnvelope);
    _supportstructures.push_back(shapeEnvelope);
  } 
  else if(geom->hasElement<mu2e::ITracker>()) 
  {
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
                                    //boost::shared_ptr<Straw> shapeDnS(new Straw(wCntPos[0],wCntPos[1],wCntPos[2], 
                                    //                                  NAN, theta, phi, wl, 
                                    //                                  _geometrymanager, _topvolume, _mainframe, 
                                    //                                  infoDnS, true));
                                    boost::shared_ptr<Straw> shapeDnS(new Straw(x,y,z, NAN, theta, phi, l, 
                                                                      _geometrymanager, _topvolume, _mainframe, 
                                                                      infoDnS, true));
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
                                    //boost::shared_ptr<Straw> shapeUpS(new Straw(wCntPos[0],wCntPos[1],wCntPos[2], 
                                    //                                  NAN, theta, phi, wl, 
                                    //                                  _geometrymanager, _topvolume, _mainframe, 
                                    //                                  infoUpS, true));
                                    boost::shared_ptr<Straw> shapeUpS(new Straw(x,y,z, NAN, theta, phi, l, 
                                                                      _geometrymanager, _topvolume, _mainframe, 
                                                                      infoUpS, true));
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
                                    boost::shared_ptr<Straw> shape(new Straw(x,y,z, NAN, theta, phi, l, 
                                                                   _geometrymanager, _topvolume, _mainframe, 
                                                                   info, true));
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
    double zHalfLength=itracker.getInnerTrackerEnvelopeParams().zHalfLength();
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
    double innerRadius=itracker.r0();//getInnerTrackerEnvelopeParams().innerRadius();
    double outerRadius=itracker.rOut();//getInnerTrackerEnvelopeParams().outerRadius();
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
    shapeEnvelope->makeGeometryVisible(true);
    _components.push_back(shapeEnvelope);
    _supportstructures.push_back(shapeEnvelope);
  }

  art::ServiceHandle<mu2e::GeometryService> geoservice;
  if(geoservice->hasElement<mu2e::DetectorSolenoid>())
  {
    mu2e::GeomHandle<mu2e::DetectorSolenoid> ds;

    double innerRadius=ds->rIn1(); 
    double outerRadius=ds->rOut2();
    double zHalfLength=ds->halfLength();
    double z=ds->position().z() - _detSysOrigin.z();

    boost::shared_ptr<ComponentInfo> infoToyDS(new ComponentInfo());
    char c[200];
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

  if(geom->hasElement<mu2e::StoppingTarget>())
  {
    mu2e::GeomHandle<mu2e::StoppingTarget> target;
    double radius=target->cylinderRadius();
    double length=target->cylinderLength();
    double z=target->centerInMu2e().z() - _detSysOrigin.z();
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
    shape->makeGeometryVisible(true);
    _components.push_back(shape);
    _supportstructures.push_back(shape);
  }

  if(geom->hasElement<mu2e::DiskCalorimeter>())
  {
    mu2e::GeomHandle<mu2e::DiskCalorimeter> calo;
    double rmax = calo->crystalHalfTrans();
    double crystalHalflength = calo->crystalHalfLength();


    int crystalIdOffset=0;
    for(unsigned int idisk=0; idisk<calo->nDisk(); idisk++)
    {
      const CLHEP::Hep3Vector &diskPos = calo->disk(idisk).origin() - _detSysOrigin;
      double innerRadius = calo->disk(idisk).size()[0];
      double outerRadius = calo->disk(idisk).size()[1];
      double diskHalflength = calo->disk(idisk).size()[2];

      findBoundaryP(_calorimeterMinmax, diskPos.x()+outerRadius, diskPos.y()+innerRadius, diskPos.z()+diskHalflength);
      findBoundaryP(_calorimeterMinmax, diskPos.x()-outerRadius, diskPos.y()-outerRadius, diskPos.z()-diskHalflength);

      char c[200];
      boost::shared_ptr<ComponentInfo> diskInfo(new ComponentInfo());
      sprintf(c,"Disk %i",idisk);
      diskInfo->setName(c);
      diskInfo->setText(0,c);
      sprintf(c,"Center at x: %.f mm, y: %.f mm, z: %.f mm",diskPos.x(),diskPos.y(),diskPos.z());
      diskInfo->setText(1,c);
      sprintf(c,"Outer radius: %.f mm, Inner radius: %.f mm, Thickness: %.f mm",outerRadius,innerRadius,2.0*diskHalflength);
      diskInfo->setText(2,c);
      boost::shared_ptr<Cylinder> calodisk(new Cylinder(diskPos.x(),diskPos.y(),diskPos.z(),  0,0,0,
                                                        diskHalflength, innerRadius, outerRadius, NAN, 
                                                        _geometrymanager, _topvolume, _mainframe, diskInfo, true));
      calodisk->makeGeometryVisible(true);
      _components.push_back(calodisk);
      _supportstructures.push_back(calodisk);

      int nCrystalInThisDisk = calo->disk(idisk).nCrystals();			
      for(int ic=0; ic<nCrystalInThisDisk; ic++)
      {
        int id=crystalIdOffset+ic;
	const CLHEP::Hep3Vector &pos = calo->crystalOrigin(crystalIdOffset+ic)-_detSysOrigin;

        boost::shared_ptr<ComponentInfo> info(new ComponentInfo());
        sprintf(c,"Disk %i, Crystal %i",idisk,id);
        info->setName(c);
        info->setText(0,c);
        sprintf(c,"Center at x: %.f mm, y: %.f mm, z: %.f mm",pos.x(),pos.y(),pos.z()+crystalHalflength);
        info->setText(1,c);
        sprintf(c,"Size: %.f mm, Thickness: %.f mm",rmax,2.0*crystalHalflength);
        info->setText(2,c);
        //these position were meant for Geant4, where the "z position of [the] hexagon is their base, not their center"
        //since this Hexagon class uses the center as a reference for, crystalHalflength needs to be added 
        boost::shared_ptr<Hexagon> shape(new Hexagon(pos.x(),pos.y(),pos.z()+crystalHalflength,
                                                     rmax,crystalHalflength,360, NAN,
                                                     _geometrymanager, _topvolume, _mainframe, info, true));
        _components.push_back(shape);
        _crystals[id]=shape;
      }
      crystalIdOffset +=nCrystalInThisDisk;
    }
  } else if(geom->hasElement<mu2e::VaneCalorimeter>())
  {
    mu2e::GeomHandle<mu2e::VaneCalorimeter> calo;
    unsigned int n=calo->nVane();
    for(unsigned int i=0; i<n; i++)
    {
      const mu2e::Vane &v=calo->vane(i);
      double x=v.origin().x() - _detSysOrigin.x();
      double y=v.origin().y() - _detSysOrigin.y();
      double z=v.origin().z() - _detSysOrigin.z();
      int    id=v.id();
      double sx=v.size().x();
      double sy=v.size().y();
      double sz=v.size().z();
      double phi=v.rotation().phi();
      double theta=v.rotation().theta();
      double psi=v.rotation().psi();

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
      boost::shared_ptr<Cube> shape(new Cube(x,y,z,  sx,sy,sz,  phi,theta,psi,   NAN,
                                        _geometrymanager, _topvolume, _mainframe, info, true));
      shape->makeGeometryVisible(true);
      _components.push_back(shape);
      _supportstructures.push_back(shape);
    }

    unsigned int roPerCrystal=calo->nROPerCrystal();
    unsigned int nro=calo->nRO();
    for(unsigned int i=0; i<nro; i+=roPerCrystal)
    {
      int vaneid=calo->vaneByRO(i);
      int crystalid=calo->crystalByRO(i);
      int rPos=calo->crystalRByRO(i);
      int zPos=calo->crystalZByRO(i);
      double crystalHalfTrans=calo->crystalHalfTrans();

      const mu2e::Vane &v=calo->vane(vaneid);
      double x=v.origin().x() - _detSysOrigin.x();
      double y=v.origin().y() - _detSysOrigin.y();
      double z=v.origin().z() - _detSysOrigin.z();
      double theta=v.rotation().theta();
      double phi=v.rotation().phi();
      double psi=v.rotation().psi();
      double sx=v.size().x();
      double sy=v.size().y();
      double sz=v.size().z();

      //Start with an unrotated vane centered at (0,0,0).
      //Before the rotation, the vector from the center of the vane
      //to the center of a crystal is (crystalX,crystalY,crystalZ).
      double crystalX=0;
      double crystalY=-sy+crystalHalfTrans*(2.0*rPos+1.0);
      double crystalZ=-sz+crystalHalfTrans*(2.0*zPos+1.0);

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
      sprintf(c,"Center at x: %.f mm, y: %.f mm, z: %.f mm",(x+rotatedX)/CLHEP::mm,(y+rotatedY)/CLHEP::mm,(z+rotatedZ)/CLHEP::mm);
      info->setText(1,c);
      sprintf(c,"Dimension  ?x: %.f mm, ?y: %.f mm, ?z: %.f mm",sx/CLHEP::mm,crystalHalfTrans/CLHEP::mm,crystalHalfTrans/CLHEP::mm);
      info->setText(2,c);
      sprintf(c,"Rotation phi: %.f °, theta: %.f °, psi: %.f °",phi/CLHEP::deg,theta/CLHEP::deg,psi/CLHEP::deg);
      info->setText(3,c);
      boost::shared_ptr<Cube> shape(new Cube(x+rotatedX,y+rotatedY,z+rotatedZ,  sx,crystalHalfTrans,crystalHalfTrans,
                                        phi,theta,psi,  NAN,
                                        _geometrymanager, _topvolume, _mainframe, info, true));
      _components.push_back(shape);
      _crystals[crystalid]=shape;
    }
  }
  //MBS
/*
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
    double mbsstartz = mbsbstsz - mbstotallen/2. - _detSysOrigin.z();
    mbsz[0] = mbsstartz + mbslen[0] ;
    mbsz[1] = mbsstartz + mbslen[0]*2. + mbslen[1];
    mbsz[2] = mbsstartz +mbstotallen - mbslen[2];
    for (unsigned int i = 0 ; i <3 ; ++i) {
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
                                               _geometrymanager, _topvolume, _mainframe, infoMBS, true));
      _components.push_back(shapeMBS);
      _mbsstructures.push_back(shapeMBS);
    }
  }
*/
  //MecoStyleProtonAbsorber 
  if(config.getBool("hasProtonAbsorber", false)) {
    if (!config.getBool("protonabsorber.isHelical", false)) {
      char c[200];
      double inr[2], outr[2], thickness, halflength, z;
      outr[0] = config.getDouble("protonabsorber.OutRadius0");
      outr[1] = config.getDouble("protonabsorber.OutRadius1");
      thickness = config.getDouble("protonabsorber.thickness");
      inr[0] = outr[0] - thickness;
      inr[1] = outr[1] - thickness;
      halflength = config.getDouble("protonabsorber.halfLength");
      mu2e::GeomHandle<mu2e::StoppingTarget> target;
      double stoppingtargetlength=target->cylinderLength();
      double stoppingtargetz=target->centerInMu2e().z() - _detSysOrigin.z();
      z = stoppingtargetz + stoppingtargetlength*0.5 + halflength;

      boost::shared_ptr<ComponentInfo> infoMecoStylePA(new ComponentInfo());
      sprintf(c,"MECOStyleProtonAbsorber");
      infoMecoStylePA->setName(c);
      infoMecoStylePA->setText(0,c);
      sprintf(c,"Inner Radius1 %.f mm  Outer Radius1 %.f mm",inr[0]/CLHEP::mm,outr[0]/CLHEP::mm);
      infoMecoStylePA->setText(1,c);
      sprintf(c,"Inner Radius2 %.f mm  Outer Radius2 %.f mm",inr[1]/CLHEP::mm,outr[1]/CLHEP::mm);
      infoMecoStylePA->setText(2,c);
      sprintf(c,"Length %.f mm",2.0*halflength/CLHEP::mm);
      infoMecoStylePA->setText(3,c);
      sprintf(c,"Center at x: 0 mm, y: 0 mm, z: %.f mm",z/CLHEP::mm);
      infoMecoStylePA->setText(4,c);
      boost::shared_ptr<Cone> shapePA(new Cone(0,0,z, 0,0,0,
                                                halflength, inr[0], outr[0], inr[1], outr[1],
                                                NAN, _geometrymanager, _topvolume, _mainframe, infoMecoStylePA, true));
      _components.push_back(shapePA);
      _mecostylepastructures.push_back(shapePA);
    }
  }

//active CRV Shields
  if( geom->hasElement<mu2e::CosmicRayShield>() ) 
  {
    mu2e::GeomHandle<mu2e::CosmicRayShield> CosmicRayShieldGeomHandle;

    std::map<std::string,mu2e::CRSScintillatorShield> const& shields = CosmicRayShieldGeomHandle->getCRSScintillatorShields();
    for(std::map<std::string,mu2e::CRSScintillatorShield>::const_iterator ishield=shields.begin(); ishield!=shields.end(); ++ishield) 
    {
      mu2e::CRSScintillatorShield const& shield = ishield->second;
      std::string shieldName = ishield->first;

      mu2e::CRSScintillatorBarDetail const& barDetail = shield.getCRSScintillatorBarDetail();
      double dx=barDetail.getHalfLengths()[0];
      double dy=barDetail.getHalfLengths()[1];
      double dz=barDetail.getHalfLengths()[2];

      int nModules = shield.nModules();
      for (int im = 0; im < nModules; ++im) 
      {
        mu2e::CRSScintillatorModule const & module = shield.getModule(im);

        int nLayers = module.nLayers();
        for (int il = 0; il < nLayers; ++il) 
        {
          mu2e::CRSScintillatorLayer const & layer = module.getLayer(il);

          int nBars = layer.nBars();
          for (int ib = 0; ib < nBars; ++ib)  
          {
            mu2e::CRSScintillatorBar const & bar = layer.getBar(ib);
            CLHEP::Hep3Vector barOffset = bar.getPosition() - _detSysOrigin;
            double x=barOffset.x();
            double y=barOffset.y();
            double z=barOffset.z();

            char c[200];
            boost::shared_ptr<ComponentInfo> info(new ComponentInfo());
            sprintf(c,"CRV Scintillator %s %i %i %i",shieldName.c_str(),im,il,ib);
            info->setName(c);
            info->setText(0,c);
            sprintf(c,"Dimension x: %.f mm, y: %.f mm, z: %.f mm",2.0*dx/CLHEP::mm,2.0*dy/CLHEP::mm,2.0*dz/CLHEP::mm);
            info->setText(1,c);
            sprintf(c,"Center at x: %.f mm, y: %.f mm, z: %.f mm",x/CLHEP::mm,y/CLHEP::mm,z/CLHEP::mm);
            info->setText(2,c);
            sprintf(c," ");
            info->setText(3,c);

            boost::shared_ptr<Cube> shape(new Cube(x,y,z,  dx,dy,dz,  0, 0, 0, NAN,
                                                   _geometrymanager, _topvolume, _mainframe, info, true));
            _components.push_back(shape);
            _crvscintillatorbars.push_back(shape);
          }
        }
      }
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
//         boost::shared_ptr<Cube> shape(new Cube(x,y,z,  dx,dy,dz,  0,0,0, NAN,
//                                           _geometrymanager, _topvolume, _mainframe, info, true));
//         _components.push_back(shape);
//         _otherstructures.push_back(shape);
//       }
//       else
//       {
// //TODO: This needs to be replaced by a shape with a hole
//         boost::shared_ptr<Cube> shape(new Cube(x,y,z,  dx,dy,dz,  0,0,0, NAN,
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
    (*structure)->makeGeometryVisible(visible);
  }

  //tracks and straws don't have to be pushed into the foreground if the structure is removed
  if(visible) toForeground();
}

void DataInterface::makeMecoStyleProtonAbsorberVisible(bool visible)
{
  std::vector<boost::shared_ptr<Cone> >::const_iterator structure;
  for(structure=_mecostylepastructures.begin(); structure!=_mecostylepastructures.end(); structure++)
  {
    (*structure)->makeGeometryVisible(visible);
  }

  //tracks and straws don't have to be pushed into the foreground if the structure is removed
  if(visible) toForeground();
}

void DataInterface::makeSupportStructuresVisible(bool visible)
{
  std::vector<boost::shared_ptr<VirtualShape> >::const_iterator structure;
  for(structure=_supportstructures.begin(); structure!=_supportstructures.end(); structure++)
  {
    (*structure)->makeGeometryVisible(visible);
  }

  //tracks and straws don't have to be pushed into the foreground if the structure is removed
  if(visible) toForeground();
}

void DataInterface::makeOtherStructuresVisible(bool visible)
{
  std::vector<boost::shared_ptr<VirtualShape> >::const_iterator structure;
  for(structure=_otherstructures.begin(); structure!=_otherstructures.end(); structure++)
  {
    (*structure)->makeGeometryVisible(visible);
  }

  //tracks and straws don't have to be pushed into the foreground if the structure is removed
  if(visible) toForeground();
}

void DataInterface::makeCrvScintillatorBarsVisible(bool visible)
{
  std::vector<boost::shared_ptr<Cube> >::const_iterator crvbars;
  for(crvbars=_crvscintillatorbars.begin(); crvbars!=_crvscintillatorbars.end(); crvbars++)
  {
    (*crvbars)->makeGeometryVisible(visible);
  }

  //tracks and straws don't have to be pushed into the foreground if the structure is removed
  if(visible) toForeground();
}

void DataInterface::toForeground()
{
  std::map<int,boost::shared_ptr<Straw> >::const_iterator straw;
  for(straw=_straws.begin(); straw!=_straws.end(); straw++)
  {
    straw->second->toForeground();
  }

  std::map<int,boost::shared_ptr<VirtualShape> >::const_iterator crystal;
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
  double mint=getHitsTimeBoundary().mint;
  double maxt=getHitsTimeBoundary().maxt;
  std::vector<boost::shared_ptr<Straw> >::const_iterator hit;
  for(hit=_hits.begin(); hit!=_hits.end(); hit++)
  {
    double time=(*hit)->getStartTime();
    if(hitcolors)
    {
      int color=TMath::FloorNint(20.0*(time-mint)/(maxt-mint));
      if(color>=20) color=19;
      if(color<=0 || isnan(color)) color=0;
      color+=2000;
      (*hit)->setColor(color);
    }
    else (*hit)->setColor(whitebackground?1:0);
  }
  std::vector<boost::shared_ptr<VirtualShape> >::const_iterator crystalhit;
  for(crystalhit=_crystalhits.begin(); crystalhit!=_crystalhits.end(); crystalhit++)
  {
    double time=(*crystalhit)->getStartTime();
    if(hitcolors)
    {
      int color=TMath::FloorNint(20.0*(time-mint)/(maxt-mint));
      if(color>=20) color=19;
      if(color<=0 || isnan(color)) color=0;
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
      int color=TMath::FloorNint(20.0*(time-mint)/(maxt-mint));
      if(color>=20) color=19;
      if(color<=0 || isnan(color)) color=0;
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
  if(steppointMChits!=nullptr)
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
  if(strawhits!=nullptr)
  {
    _numberHits=strawhits->size();
    std::vector<mu2e::StrawHit>::const_iterator iter;
    int hitnumber=0;
    for(iter=strawhits->begin(); iter!=strawhits->end(); iter++, hitnumber++)
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
          straw->second->setHitNumber(hitnumber);
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
  const mu2e::KalRepCollection *kalRepHits=contentSelector->getSelectedHitCollection<mu2e::KalRepCollection>();
  if(kalRepHits!=nullptr)
  {
    boost::shared_ptr<TGraphErrors> residualGraph(new TGraphErrors());
    residualGraph->SetTitle("Residual Graph");

    for(unsigned int i=0; i<kalRepHits->size(); i++)
    {
      const KalRep &particle = kalRepHits->at(i);
      const TrkHotList* hots = particle.hotList();
      if(hots!=nullptr)
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
            double hitT0 = strawHit->hitT0()._t0; //this is the time the hit "arrived at the straw"
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

  const mu2e::StepPointMCCollection *calosteppoints=contentSelector->getSelectedCaloHitCollection<mu2e::StepPointMCCollection>();
  if(calosteppoints!=nullptr)
  {
    _numberCrystalHits=calosteppoints->size();
    std::vector<mu2e::StepPointMC>::const_iterator iter;
    for(iter=calosteppoints->begin(); iter!=calosteppoints->end(); iter++)
    {
      const mu2e::StepPointMC& calohit = *iter;
      int crystalid = calohit.volumeId();
      int trackid = calohit.trackId().asInt();
      double time = calohit.time();
      double energy = calohit.eDep();
      std::map<int,boost::shared_ptr<VirtualShape> >::iterator crystal=_crystals.find(crystalid);
      if(crystal!=_crystals.end() && !isnan(time))
      {
        double previousStartTime=crystal->second->getStartTime();
        if(isnan(previousStartTime))
        {
          findBoundaryT(_hitsTimeMinmax, time);  //is it Ok to exclude all following hits from the time window?
          crystal->second->setStartTime(time);
          char c[100];
          sprintf(c,"hit time(s): %gns",time/CLHEP::ns);
          crystal->second->getComponentInfo()->setText(2,c);
          sprintf(c,"deposited energy(s): %geV",energy/CLHEP::eV);
          crystal->second->getComponentInfo()->setText(3,c);
          sprintf(c,"track ID(s): %i",trackid);
          crystal->second->getComponentInfo()->setText(4,c);
          _crystalhits.push_back(crystal->second);
        }
        else
        {
          crystal->second->getComponentInfo()->expandLine(2,time/CLHEP::ns,"ns");
          crystal->second->getComponentInfo()->expandLine(3,energy/CLHEP::eV,"eV");
          crystal->second->getComponentInfo()->expandLine(4,trackid,"");
        }
      }
    }
  }

  const mu2e::CaloCrystalHitCollection *calocrystalhits=contentSelector->getSelectedCaloHitCollection<mu2e::CaloCrystalHitCollection>();
  if(calocrystalhits!=nullptr)
  {
    _numberCrystalHits=calocrystalhits->size();
    std::vector<mu2e::CaloCrystalHit>::const_iterator iter;
    for(iter=calocrystalhits->begin(); iter!=calocrystalhits->end(); iter++)
    {
      const mu2e::CaloCrystalHit& calohit = *iter;
      int crystalid = calohit.id();
      double time = calohit.time();
      double energy = calohit.energyDep();
      std::map<int,boost::shared_ptr<VirtualShape> >::iterator crystal=_crystals.find(crystalid);
      if(crystal!=_crystals.end() && !isnan(time))
      {
        double previousStartTime=crystal->second->getStartTime();
        if(isnan(previousStartTime))
        {
          findBoundaryT(_hitsTimeMinmax, time);  //is it Ok to exclude all following hits from the time window?
          crystal->second->setStartTime(time);
          char c[100];
          sprintf(c,"hit time(s): %gns",time/CLHEP::ns);
          crystal->second->getComponentInfo()->setText(2,c);
          sprintf(c,"deposited energy(s): %geV",energy/CLHEP::eV);
          crystal->second->getComponentInfo()->setText(3,c);
          _crystalhits.push_back(crystal->second);
        }
        else
        {
          crystal->second->getComponentInfo()->expandLine(2,time/CLHEP::ns,"ns");
          crystal->second->getComponentInfo()->expandLine(3,energy/CLHEP::eV,"eV");
        }
      }
    }
  }

  const mu2e::CaloHitCollection *calohits=contentSelector->getSelectedCaloHitCollection<mu2e::CaloHitCollection>();
  art::ServiceHandle<mu2e::GeometryService> geoservice;
  if(calohits!=nullptr && geoservice->hasElement<mu2e::VaneCalorimeter>())
  {
    mu2e::GeomHandle<mu2e::VaneCalorimeter> calo;
    _numberCrystalHits=calohits->size();  //this is not accurate since the return value gives the RO hits
    std::vector<mu2e::CaloHit>::const_iterator iter;
    for(iter=calohits->begin(); iter!=calohits->end(); iter++)
    {
      const mu2e::CaloHit& calohit = *iter;
      int roid = calohit.id();
      int crystalid=calo->crystalByRO(roid);
      double time = calohit.time();
      double energy = calohit.energyDep();
      std::map<int,boost::shared_ptr<VirtualShape> >::iterator crystal=_crystals.find(crystalid);
      if(crystal!=_crystals.end() && !isnan(time))
      {
        double previousStartTime=crystal->second->getStartTime();
        if(isnan(previousStartTime))
        {
          findBoundaryT(_hitsTimeMinmax, time);  //is it Ok to exclude all following hits from the time window?
          crystal->second->setStartTime(time);
          char c[100];
          sprintf(c,"hit time(s): %gns",time/CLHEP::ns);
          crystal->second->getComponentInfo()->setText(2,c);
          sprintf(c,"deposited energy(s): %geV",energy/CLHEP::eV);
          crystal->second->getComponentInfo()->setText(3,c);
          sprintf(c,"RO ID(s): %i",roid);
          crystal->second->getComponentInfo()->setText(4,c);
          _crystalhits.push_back(crystal->second);
        }
        else
        {
          crystal->second->getComponentInfo()->expandLine(2,time/CLHEP::ns,"ns");
          crystal->second->getComponentInfo()->expandLine(3,energy/CLHEP::eV,"eV");
          crystal->second->getComponentInfo()->expandLine(4,roid,"");
        }
      }
    }
  }

  unsigned int physicalVolumeEntries=0;
  const mu2e::PhysicalVolumeInfoCollection *physicalVolumes=contentSelector->getPhysicalVolumeInfoCollection();
  const mu2e::PhysicalVolumeInfoMultiCollection *physicalVolumesMulti=contentSelector->getPhysicalVolumeInfoMultiCollection();
  if(physicalVolumes!=nullptr)
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
      int id = particle.id().asInt();   //is identical with cet::map_vector_key& particleKey = iter->first;
      int parentid = particle.parentId().asInt();
      int particleid=particle.pdgId();
      int trackclass=trackInfos[i].classID;
      int trackclassindex=trackInfos[i].index;
      std::string particlecollection=trackInfos[i].entryText;
      std::string particlename=HepPID::particleName(particle.pdgId());
      std::string startVolumeName="unknown volume";
      std::string endVolumeName="unknown volume";
      if(physicalVolumes!=nullptr)
      {
        unsigned int startVolume=particle.startVolumeIndex();
        unsigned int endVolume  =particle.endVolumeIndex();
        if(startVolume<physicalVolumeEntries && startVolume>=0) startVolumeName=physicalVolumes->at(startVolume).name();
        if(endVolume<physicalVolumeEntries && endVolume>=0) endVolumeName=physicalVolumes->at(endVolume).name();
      }
      else if(physicalVolumesMulti!=nullptr)
      {
        mu2e::PhysicalVolumeMultiHelper volumeMultiHelper(*physicalVolumesMulti);
        startVolumeName=volumeMultiHelper.startVolume(particle).name();
        endVolumeName=volumeMultiHelper.endVolume(particle).name();
      }
      double x1=particle.startPosition().x() - _detSysOrigin.x();
      double y1=particle.startPosition().y() - _detSysOrigin.y();
      double z1=particle.startPosition().z() - _detSysOrigin.z();
      double t1=particle.startGlobalTime();
      double e1=particle.startMomentum().e();
      double x2=particle.endPosition().x() - _detSysOrigin.x();
      double y2=particle.endPosition().y() - _detSysOrigin.y();
      double z2=particle.endPosition().z() - _detSysOrigin.z();
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
      boost::shared_ptr<Track> shape(new Track(x1,y1,z1,t1, x2,y2,z2,t2, 
                                               particleid, trackclass, trackclassindex, e1, 
                                               _geometrymanager, _topvolume, _mainframe, info, false));
      findTrajectory(contentSelector,shape,particle.id(), t1,t2, simParticles,particle.daughterIds(), trackInfos[i]);
      _components.push_back(shape);
      _tracks.push_back(shape);
    }
  }

#ifdef BABARINSTALLED
  trackInfos.clear();
  std::vector<const mu2e::KalRepCollection*> kalRepCollectionVector=contentSelector->getSelectedTrackCollection<mu2e::KalRepCollection>(trackInfos);
  for(unsigned int i=0; i<kalRepCollectionVector.size(); i++)
  {
    const mu2e::KalRepCollection *kalReps=kalRepCollectionVector[i];
    for(unsigned int j=0; j<kalReps->size(); j++)
    {
      KalRep const* kalrep = kalReps->get(j);
      double t0=kalrep->t0().t0();
      {
        int trackclass=trackInfos[i].classID;
        int trackclassindex=trackInfos[i].index;
        std::string trackcollection=trackInfos[i].entryText;
        int particleid=kalrep->particleType().particleType();
        std::string particlename=HepPID::particleName(particleid);
        char c0[200], c2[200], c3[200], c4[200];
        sprintf(c0,"Kalman Track %i  %s  (%s)",j,particlename.c_str(),trackcollection.c_str());
        boost::shared_ptr<ComponentInfo> info(new ComponentInfo());
        info->setName(c0);
        info->setText(0,c0);

        double hitcount=0;
        double offset=0;
        const TrkHotList* hots=kalrep->hotList();
        if(hots!=nullptr)
        {
          boost::shared_ptr<TGraphErrors> residualGraph(new TGraphErrors());
          residualGraph->SetTitle("Residual Graph");
          for(TrkHotList::hot_iterator iter=hots->begin(); iter!=hots->end(); iter++)
          {
            const TrkHitOnTrk *hitOnTrack = iter.get();
            const mu2e::TrkStrawHit* strawHit = dynamic_cast<const mu2e::TrkStrawHit*>(hitOnTrack);
            if(strawHit)
            {
              double strawTime   = strawHit->hitT0()._t0/CLHEP::ns;
              double trackTime   = strawTime*CLHEP::ns;  //TODO: add correction for drift time
              double weight= strawHit->weight();  
              double fltLen= strawHit->fltLen();
              const HepPoint &p=strawHit->hitTraj()->position(strawHit->hitLen());
              double t     = kalrep->arrivalTime(fltLen);
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

        double fltLMin=kalrep->startValidRange();
        double fltLMax=kalrep->endValidRange();
        double p1=kalrep->momentum(fltLMin).mag();
        double p2=kalrep->momentum(fltLMax).mag();
        double x1=kalrep->position(fltLMin).x();
        double y1=kalrep->position(fltLMin).y();
        double z1=kalrep->position(fltLMin).z();
        double x2=kalrep->position(fltLMax).x();
        double y2=kalrep->position(fltLMax).y();
        double z2=kalrep->position(fltLMax).z();
        double t1=kalrep->arrivalTime(fltLMin)+offset;
        double t2=kalrep->arrivalTime(fltLMax)+offset;
        boost::shared_ptr<Track> track(new Track(x1,y1,z1,t1, x2,y2,z2,t2, 
                                                 particleid, trackclass, trackclassindex, p1, 
                                                 _geometrymanager, _topvolume, _mainframe, info, false));
        _components.push_back(track);
        _tracks.push_back(track);

        double fltStep = (fltLMax - fltLMin)/400.0;
        for(unsigned int step = 0; step <= 400.0; step++) 
        {		
          double fltL = fltLMin + step*fltStep;
          double   t = kalrep->arrivalTime(fltL)+offset;
          HepPoint p = kalrep->position(fltL);
          findBoundaryT(_tracksTimeMinmax, t);
          findBoundaryP(_tracksMinmax, p.x(), p.y(), p.z());
          track->addTrajectoryPoint(p.x(), p.y(), p.z(), t);
        }
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

  // TrkExt track
  trackInfos.clear();
  std::vector<const mu2e::TrkExtTrajCollection*> trkExtTrajCollectionVector=contentSelector->getSelectedTrackCollection<mu2e::TrkExtTrajCollection>(trackInfos);
  for(unsigned int i=0; i<trkExtTrajCollectionVector.size(); i++)
  {
    // Read a TrkExtTrajCollection
    const mu2e::TrkExtTrajCollection & trkExtTrajCollection = *trkExtTrajCollectionVector[i];
    for(unsigned int j=0; j<trkExtTrajCollection.size(); j++)
    {
      // read a TrkExtTraj
      const mu2e::TrkExtTraj &trkExtTraj = trkExtTrajCollection.at(j);
      int particleid=11;
      int trackclass=trackInfos[i].classID;
      int trackclassindex=trackInfos[i].index;
      std::string trackcollection=trackInfos[i].entryText;

      std::string particlename=HepPID::particleName(particleid);
      char c0[200];
      sprintf(c0,"TrkExt Trajectory %i  %s  (%s)", trkExtTraj.id(), particlename.c_str(),trackcollection.c_str());
      boost::shared_ptr<ComponentInfo> info(new ComponentInfo());
      info->setName(c0);
      info->setText(0,c0);

      double p1 = trkExtTraj.front().momentum().mag();
      double x1 = trkExtTraj.front().x();
      double y1 = trkExtTraj.front().y();
      double z1 = trkExtTraj.front().z();
      double x2 = trkExtTraj.back().x();
      double y2 = trkExtTraj.back().y();
      double z2 = trkExtTraj.back().z();
      double t1 = 0;
      double t2 = 0;
      boost::shared_ptr<Track> track(new Track(x1,y1,z1,t1, x2,y2,z2,t2, 
                                               particleid, trackclass, trackclassindex, p1, 
                                               _geometrymanager, _topvolume, _mainframe, info, false));
      _components.push_back(track);
      _tracks.push_back(track);

      for (unsigned int k = 0 ; k < trkExtTraj.size() ; k+=10) {
        const mu2e::TrkExtTrajPoint & trkExtTrajPoint = trkExtTraj[k];
        track->addTrajectoryPoint(trkExtTrajPoint.x(), trkExtTrajPoint.y(), trkExtTrajPoint.z(), 0);
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
  const mu2e::MCTrajectoryCollection *mcTrajectories=contentSelector->getMCTrajectoryCollection(trackInfo);
  if(mcTrajectories!=nullptr)
  {
    std::map<art::Ptr<mu2e::SimParticle>,mu2e::MCTrajectory>::const_iterator traj_iter;
    for(traj_iter=mcTrajectories->begin(); traj_iter!=mcTrajectories->end(); traj_iter++)
    {
      if(traj_iter->first->id()==id) 
      {
        const std::vector<CLHEP::HepLorentzVector> &points = traj_iter->second.points();
        std::vector<CLHEP::HepLorentzVector>::const_iterator point_iter;
        for(point_iter=points.begin(); point_iter!=points.end(); point_iter++)
        {
          track->addTrajectoryPoint(point_iter->x()-_detSysOrigin.x(),
                                    point_iter->y()-_detSysOrigin.y(),
                                    point_iter->z()-_detSysOrigin.z(),
                                    point_iter->t());
        }
      }
    }
    return;  //otherwise try point trajectories
  }


#ifdef USETRAJECTORY
  const mu2e::PointTrajectoryCollection *pointTrajectories=contentSelector->getPointTrajectoryCollection(trackInfo);
  if(pointTrajectories!=nullptr)
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
          ts.v = pVect[i] - _detSysOrigin;
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
              CLHEP::Hep3Vector v(particle->startPosition() - _detSysOrigin);
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
    if(!(*iter)->isGeometry()) {iter=_components.erase(iter);} //things like tracks and drift radii
    else iter++;
  }

  std::vector<boost::shared_ptr<Straw> >::const_iterator hit;
  for(hit=_hits.begin(); hit!=_hits.end(); hit++)
  {
    for(int i=1; i<5; i++) (*hit)->getComponentInfo()->removeLine(i);  //keep first line
    (*hit)->getComponentInfo()->getHistVector().clear();
    (*hit)->setHitNumber(-1);
    (*hit)->setStartTime(NAN);
    (*hit)->start();
  }
  std::vector<boost::shared_ptr<VirtualShape> >::const_iterator crystalhit;
  for(crystalhit=_crystalhits.begin(); crystalhit!=_crystalhits.end(); crystalhit++)
  {
    for(int i=1; i<5; i++) (*crystalhit)->getComponentInfo()->removeLine(i); //keep first line
    (*crystalhit)->getComponentInfo()->getHistVector().clear();
    (*crystalhit)->setStartTime(NAN);
    (*crystalhit)->start();
  }

  _hits.clear();
  _crystalhits.clear();
  _tracks.clear();  //will call the d'tors of all tracks, since they aren't used anywhere anymore
  _driftradii.clear(); //will call the d'tors of all driftradii, since they aren't used anywhere anymore

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
  _crvscintillatorbars.clear();
  _mbsstructures.clear();
  _mecostylepastructures.clear();
  delete _geometrymanager;
  _geometrymanager=nullptr;

  _mainframe->getHistDrawVector().clear();
}

}
