//
using namespace std;
#include "Offline/EventDisplay/src/DataInterface.h"

#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/CRVConditions/inc/CRVDigitizationPeriod.hh"
#include "Offline/DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "Offline/EventDisplay/src/Cube.h"
#include "Offline/EventDisplay/src/Cylinder.h"
#include "Offline/EventDisplay/src/Cone.h"
#include "Offline/EventDisplay/src/EventDisplayFrame.h"
#include "Offline/EventDisplay/src/Hexagon.h"
#include "Offline/EventDisplay/src/Straw.h"
#include "Offline/EventDisplay/src/Track.h"
#include "Offline/EventDisplay/src/TrackColorSelector.h"
#include "Offline/EventDisplay/src/dict_classes/ComponentInfo.h"
#include "Offline/EventDisplay/src/dict_classes/EventDisplayViewSetup.h"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/DataProducts/inc/GenVector.hh"
#include "Offline/MCDataProducts/inc/PhysicalVolumeInfoMultiCollection.hh"
#include "Offline/MCDataProducts/inc/MCTrajectoryCollection.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/Mu2eUtilities/inc/PhysicalVolumeMultiHelper.hh"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/CrvDigi.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/TrkExtTraj.hh"
#include "Offline/StoppingTargetGeom/inc/StoppingTarget.hh"
#include "Offline/StoppingTargetGeom/inc/TargetFoil.hh"
#include "Offline/TrkReco/inc/TrkUtilities.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "art/Framework/Principal/Run.h"
#include "cetlib/map_vector.h"
#include <TAxis3D.h>
#include <TGFrame.h>
#include <TGeoVolume.h>
#include <TMath.h>
#include <TView.h>
#include <TGraphErrors.h>
#include <TF1.h>

#include <boost/shared_array.hpp>

using namespace CLHEP;
//#include "Offline/RecoDataProducts/inc/KalRepCollection.hh"
//#include "Offline/BTrkData/inc/TrkStrawHit.hh"
//#include "BTrk/KalmanTrack/KalRep.hh"
#include "Offline/BTrkLegacy/inc/ExternalInfo.hh"

namespace mu2e_eventdisplay
{

DataInterface::DataInterface(EventDisplayFrame *mainframe, double kalStepSize):
              _geometrymanager(nullptr),_topvolume(nullptr),_mainframe(mainframe),_kalStepSize(kalStepSize)
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

    _particleInfo = make_unique<mu2e::ParticleInfo>();
    ExternalInfo::set(_particleInfo.get());
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

  std::vector<boost::shared_ptr<Cube> >::const_iterator crvhit;
  for(crvhit=_crvhits.begin(); crvhit!=_crvhits.end(); crvhit++)
  {
    (*crvhit)->setFilter(_minTime, _maxTime);
    (*crvhit)->update(time);
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
   TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0,0,0);
   TGeoMedium *Vacuum = new TGeoMedium("Vacuum",1, matVacuum);
  _topvolume = _geometrymanager->MakeBox("TopVolume", Vacuum, 1000, 1000, 1500);
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
  resetBoundaryP(_crvMinmax);

  art::ServiceHandle<mu2e::GeometryService> geom;

  const mu2e::SimpleConfig &config = geom->config();
  _detSysOrigin = mu2e::GeomHandle<mu2e::DetectorSystem>()->getOrigin();

  if(geom->hasElement<mu2e::Tracker>())
  {
    mu2e::GeomHandle<mu2e::Tracker> tracker;
    double trackerOffset = tracker->g4Tracker()->z0() - _detSysOrigin.z();

//Straws
    const auto& allStraws = tracker->getStraws();
    // for(const auto & elem : allStraws)
    for (size_t i = 0; i<tracker->nStraws(); ++i)
    {
      // const mu2e::Straw& s = elem;
      const mu2e::Straw& s = allStraws[i];
      const CLHEP::Hep3Vector& p = s.getMidPoint();
      const CLHEP::Hep3Vector& d = s.getDirection();
      double x = p.x();
      double y = p.y();
      double z = p.z()+trackerOffset;
      double theta = d.theta();
      double phi = d.phi();
      double l = s.halfLength();
      int idStraw =  s.id().getStraw();
      int idLayer =  s.id().getLayer();
      int idPanel =  s.id().getPanel();
      int idPlane =  s.id().getPlane();
      int id = s.id().asUint16();

      boost::shared_ptr<ComponentInfo> info(new ComponentInfo());
      std::string c=Form("Straw %i  Layer %i  Panel %i  Plane %i",idStraw,idLayer,idPanel,idPlane);
      info->setName(c.c_str());
      info->setText(0,c.c_str());
      boost::shared_ptr<Straw> shape(new Straw(x,y,z, NAN, theta, phi, l,
                                               _geometrymanager, _topvolume, _mainframe, info, true));
      _components.push_back(shape);
      _straws[id]=shape;
    }

//Support Structure
    double innerRadius=tracker->g4Tracker()->getSupportParams().innerRadius();
    double outerRadius=tracker->g4Tracker()->getSupportParams().outerRadius();
    double zHalfLength=tracker->g4Tracker()->getInnerTrackerEnvelopeParams().zHalfLength();
    findBoundaryP(_trackerMinmax, outerRadius, outerRadius, zHalfLength);
    findBoundaryP(_trackerMinmax, -outerRadius, -outerRadius, -zHalfLength);

    boost::shared_ptr<ComponentInfo> info(new ComponentInfo());
    info->setName("Tracker Support Structure");
    info->setText(0,"Tracker Support Structure");
    info->setText(1,Form("Inner Radius %.f mm  Outer Radius %.f mm",innerRadius/CLHEP::mm,outerRadius/CLHEP::mm));
    info->setText(2,Form("Length %.f mm",2.0*zHalfLength/CLHEP::mm));
    info->setText(3,Form("Center at x: 0 mm, y: 0 mm, z: %.f mm",trackerOffset));
    boost::shared_ptr<Cylinder> shape(new Cylinder(0,0,trackerOffset, 0,0,trackerOffset,
                                          zHalfLength,innerRadius,outerRadius, NAN,
                                          _geometrymanager, _topvolume, _mainframe, info, true));
    shape->makeGeometryVisible(true);
    _components.push_back(shape);
    _supportstructures.push_back(shape);

//Envelope
    innerRadius=tracker->g4Tracker()->getInnerTrackerEnvelopeParams().innerRadius();
    outerRadius=tracker->g4Tracker()->getInnerTrackerEnvelopeParams().outerRadius();
    zHalfLength=tracker->g4Tracker()->getInnerTrackerEnvelopeParams().zHalfLength();

    boost::shared_ptr<ComponentInfo> infoEnvelope(new ComponentInfo());
    infoEnvelope->setName("Tracker Envelope");
    infoEnvelope->setText(0,"Tracker Envelope");
    infoEnvelope->setText(1,Form("Inner Radius %.f mm  Outer Radius %.f mm",innerRadius/CLHEP::mm,outerRadius/CLHEP::mm));
    infoEnvelope->setText(2,Form("Length %.f mm",2.0*zHalfLength/CLHEP::mm));
    infoEnvelope->setText(3,Form("Center at x: 0 mm, y: 0 mm, z: %.f mm",trackerOffset));
    boost::shared_ptr<Cylinder> shapeEnvelope(new Cylinder(0,0,trackerOffset, 0,0,trackerOffset,
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
    infoToyDS->setName("Toy DS");
    infoToyDS->setText(0,"Toy DS");
    infoToyDS->setText(1,Form("Inner Radius %.f mm  Outer Radius %.f mm",innerRadius/CLHEP::mm,outerRadius/CLHEP::mm));
    infoToyDS->setText(2,Form("Length %.f mm",2.0*zHalfLength/CLHEP::mm));
    infoToyDS->setText(3,Form("Center at x: 0 mm, y: 0 mm, z: %.f mm",z/CLHEP::mm));
    boost::shared_ptr<Cylinder> shapeToyDS(new Cylinder(0,0,z, 0,0,0,
                                               zHalfLength,innerRadius,outerRadius, NAN,
                                               _geometrymanager, _topvolume, _mainframe, infoToyDS, true));
    _components.push_back(shapeToyDS);
    _otherstructures.push_back(shapeToyDS);
  }

  if(geom->hasElement<mu2e::StoppingTarget>())
  {
    mu2e::GeomHandle<mu2e::StoppingTarget> target;
    unsigned int n=target->nFoils();
    for(unsigned int i=0; i<n; i++)
    {
      const mu2e::TargetFoil &foil=target->foil(i);
      int id = foil.id();
      double x = foil.centerInDetectorSystem().x();
      double y = foil.centerInDetectorSystem().y();
      double z = foil.centerInDetectorSystem().z();
      double radius = foil.rOut();
      double halfThickness = foil.halfThickness();

      findBoundaryP(_targetMinmax, x+radius, y+radius, z+halfThickness);
      findBoundaryP(_targetMinmax, x-radius, y-radius, z-halfThickness);

      boost::shared_ptr<ComponentInfo> info(new ComponentInfo());
      std::string c=Form("Target Foil ID %i",id);
      info->setName(c.c_str());
      info->setText(0,c.c_str());
      info->setText(1,Form("Outer Radius %.2f mm",radius/CLHEP::mm));
      info->setText(2,Form("Half thickness %.2f mm",halfThickness/CLHEP::mm));
      info->setText(3,Form("Center at x: %.2f mm, y: %.2f mm, z: %.2f mm",x/CLHEP::mm,y/CLHEP::mm,z/CLHEP::mm));
      boost::shared_ptr<Cylinder> shape(new Cylinder(x,y,z, 0,0,0, halfThickness,0,radius, NAN,
                                          _geometrymanager, _topvolume, _mainframe, info, true));
      shape->makeGeometryVisible(true);
      _components.push_back(shape);
      _supportstructures.push_back(shape);
    }
  }

  if(geom->hasElement<mu2e::DiskCalorimeter>())
  {
    mu2e::GeomHandle<mu2e::DiskCalorimeter> calo;

    double diskCaseDZLength      = calo->caloInfo().getDouble("diskCaseZLength")/2.0;
    double diskRadiusIn          = calo->caloInfo().getDouble("caloDiskRadiusIn");
    double diskRadiusOut         = calo->caloInfo().getDouble("caloDiskRadiusOut");

    double FPCarbonDZ               = calo->caloInfo().getDouble("FPCarbonZLength")/2.0;
    double FPFoamDZ                 = calo->caloInfo().getDouble("FPFoamZLength")/2.0;
    double FPCoolPipeRadius         = calo->caloInfo().getDouble("FPCoolPipeRadius");
    double pipeRadius               = calo->caloInfo().getDouble("pipeRadius");
    double frontPanelHalfThick      = (2.0*FPCarbonDZ+2.0*FPFoamDZ-pipeRadius+FPCoolPipeRadius)/2.0;

    double crystalDXY            = calo->caloInfo().getDouble("crystalXYLength")/2.0;
    double crystalDZ             = calo->caloInfo().getDouble("crystalZLength")/2.0;
    double crystalFrameDZ        = calo->caloInfo().getDouble("crystalCapZLength")/2.0;
    double wrapperHalfThick      = calo->caloInfo().getDouble("wrapperThickness")/2.0;
    double wrapperDXY            = crystalDXY + 2.0*wrapperHalfThick;
    double wrapperDZ             = crystalDZ + 2.0*crystalFrameDZ;

    double FEEDZ                = calo->caloInfo().getDouble("FEEZLength")/2.0;
    double FEEBoxThickness      = calo->caloInfo().getDouble("FEEBoxThickness");
    double FEEBoxDZ             = FEEDZ + 2*FEEBoxThickness;
    double BPPipeRadiusHigh     = calo->caloInfo().getDouble("BPPipeRadiusHigh");
    double BPPipeDZOffset       = calo->caloInfo().getDouble("BPPipeZOffset")/2.0;
    double BPFEEDZ              = FEEBoxDZ + BPPipeDZOffset + BPPipeRadiusHigh;
    double holeDZ               = calo->caloInfo().getDouble("BPHoleZLength")/2.0;
    double zHalfBP              = BPFEEDZ+holeDZ;

    double crystalDiskLogOffset = frontPanelHalfThick - zHalfBP;

    int icrystal=0;
    for(unsigned int idisk=0; idisk<calo->nDisks(); idisk++)
    {
      CLHEP::Hep3Vector diskPos = calo->disk(idisk).geomInfo().origin() - _detSysOrigin;
      diskPos += CLHEP::Hep3Vector(0.0, 0.0, crystalDiskLogOffset);

      findBoundaryP(_calorimeterMinmax, diskPos.x()+diskRadiusOut, diskPos.y()+diskRadiusOut, diskPos.z()+diskCaseDZLength);
      findBoundaryP(_calorimeterMinmax, diskPos.x()-diskRadiusOut, diskPos.y()-diskRadiusOut, diskPos.z()-diskCaseDZLength);

      boost::shared_ptr<ComponentInfo> diskInfo(new ComponentInfo());
      std::string c=Form("Disk %i",idisk);
      diskInfo->setName(c.c_str());
      diskInfo->setText(0,c.c_str());
      diskInfo->setText(1,Form("Center at x: %.f mm, y: %.f mm, z: %.f mm",diskPos.x(),diskPos.y(),diskPos.z()));
      diskInfo->setText(2,Form("Outer radius: %.f mm, Inner radius: %.f mm, Thickness: %.f mm",diskRadiusOut,diskRadiusIn,2.0*diskCaseDZLength));
      boost::shared_ptr<Cylinder> calodisk(new Cylinder(diskPos.x(),diskPos.y(),diskPos.z(),  0,0,0,
                                                        diskCaseDZLength, diskRadiusIn, diskRadiusOut, NAN,
                                                        _geometrymanager, _topvolume, _mainframe, diskInfo, true));
      calodisk->makeGeometryVisible(true);
      _components.push_back(calodisk);
      _supportstructures.push_back(calodisk);

      int nCrystalInThisDisk = calo->disk(idisk).nCrystals();
      for(int ic=0; ic<nCrystalInThisDisk; ic++)
      {
        CLHEP::Hep3Vector crystalPosition = calo->disk(idisk).crystal(ic).localPosition();
        crystalPosition.setZ(diskCaseDZLength-wrapperDZ);
        crystalPosition += diskPos;

        boost::shared_ptr<ComponentInfo> info(new ComponentInfo());
        std::string c=Form("Disk %i, Crystal %i",idisk,ic);
        info->setName(c.c_str());
        info->setText(0,c.c_str());
        info->setText(1,Form("Center at x: %.f mm, y: %.f mm, z: %.f mm",crystalPosition.x(),crystalPosition.y(),crystalPosition.z()));
        info->setText(2,Form("XYSize: %.f mm, Thickness: %.f mm",2.0*wrapperDXY,2.0*wrapperDZ));
        boost::shared_ptr<Hexagon> shape(new Hexagon(crystalPosition.x(),crystalPosition.y(),crystalPosition.z(),
                                                     wrapperDXY,wrapperDZ,360, NAN,
                                                     _geometrymanager, _topvolume, _mainframe, info, true));
        _components.push_back(shape);
        _crystals[icrystal]=shape;
        icrystal++;
      }
    }
  }

  //MBS
/*
  if(config.getBool("hasMBS", false))
  {
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
    for (unsigned int i = 0 ; i <3 ; ++i)
    {
      boost::shared_ptr<ComponentInfo> infoMBS(new ComponentInfo());
      std::string c=Form(c,"MBS %d", i);
      infoMBS->setName(c.c_str());
      infoMBS->setText(0,c.c_str());
      infoMBS->setText(1,Form("Inner Radius %.f mm  Outer Radius %.f mm",mbsinr[i]/CLHEP::mm,mbsoutr[i]/CLHEP::mm));
      infoMBS->setText(2,Form("Length %.f mm",2.0*mbslen[i]/CLHEP::mm));
      infoMBS->setText(3,Form("Center at x: 0 mm, y: 0 mm, z: %.f mm",mbsz[i]/CLHEP::mm));
      boost::shared_ptr<Cylinder> shapeMBS(new Cylinder(0,0,mbsz[i], 0,0,0,
                                               mbslen[i], mbsinr[i], mbsoutr[i], NAN,
                                               _geometrymanager, _topvolume, _mainframe, infoMBS, true));
      _components.push_back(shapeMBS);
      _mbsstructures.push_back(shapeMBS);
    }
  }
*/
  //MecoStyleProtonAbsorber
  if(config.getBool("hasProtonAbsorber", false))
  {
    if (!config.getBool("protonabsorber.isHelical", false))
    {
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
      infoMecoStylePA->setName("MECOStyleProtonAbsorber");
      infoMecoStylePA->setText(0,"MECOStyleProtonAbsorber");
      infoMecoStylePA->setText(1,Form("Inner Radius1 %.f mm  Outer Radius1 %.f mm",inr[0]/CLHEP::mm,outr[0]/CLHEP::mm));
      infoMecoStylePA->setText(2,Form("Inner Radius2 %.f mm  Outer Radius2 %.f mm",inr[1]/CLHEP::mm,outr[1]/CLHEP::mm));
      infoMecoStylePA->setText(3,Form("Length %.f mm",2.0*halflength/CLHEP::mm));
      infoMecoStylePA->setText(4,Form("Center at x: 0 mm, y: 0 mm, z: %.f mm",z/CLHEP::mm));
      boost::shared_ptr<Cone> shapePA(new Cone(0,0,z, 0,0,0,
                                                halflength, inr[0], outr[0], inr[1], outr[1],
                                                NAN, _geometrymanager, _topvolume, _mainframe, infoMecoStylePA, true));
      _components.push_back(shapePA);
      _mecostylepastructures.push_back(shapePA);
    }
  }

//CRV
  if( geom->hasElement<mu2e::CosmicRayShield>() )
  {
    mu2e::GeomHandle<mu2e::CosmicRayShield> CosmicRayShieldGeomHandle;
    std::vector<mu2e::CRSScintillatorShield> const& shields = CosmicRayShieldGeomHandle->getCRSScintillatorShields();
    for(std::vector<mu2e::CRSScintillatorShield>::const_iterator ishield=shields.begin(); ishield!=shields.end(); ++ishield)
    {
      mu2e::CRSScintillatorShield const& shield = *ishield;
      std::string const& shieldName = shield.getName();
//if(shieldName.find("CRV_C1")!=std::string::npos) continue;
//if(shieldName.find("CRV_C2")!=std::string::npos) continue;
//if(shieldName.find("CRV_T")!=std::string::npos) continue;

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
            int index = bar.index().asInt();
            CLHEP::Hep3Vector barOffset = bar.getPosition() - _detSysOrigin;
            double x=barOffset.x();
            double y=barOffset.y();
            double z=barOffset.z();

            findBoundaryP(_crvMinmax, x+dx, y+dy, z+dz);
            findBoundaryP(_crvMinmax, x-dx, y-dy, z-dz);

            boost::shared_ptr<ComponentInfo> info(new ComponentInfo());
            std::string c=Form("CRV Scintillator %s  module %i  layer %i  bar %i  (index %i)",shieldName.c_str(),im,il,ib, index);
            info->setName(c.c_str());
            info->setText(0,c.c_str());
            info->setText(1,Form("Dimension x: %.f mm, y: %.f mm, z: %.f mm",2.0*dx/CLHEP::mm,2.0*dy/CLHEP::mm,2.0*dz/CLHEP::mm));
            info->setText(2,Form("Center at x: %.f mm, y: %.f mm, z: %.f mm",x/CLHEP::mm,y/CLHEP::mm,z/CLHEP::mm));

            boost::shared_ptr<Cube> shape(new Cube(x,y,z,  dx,dy,dz,  0, 0, 0, NAN,
                                                   _geometrymanager, _topvolume, _mainframe, info, true));
            _components.push_back(shape);
            _crvscintillatorbars[index]=(shape);
          }
        }
      }
    }
  }
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
  std::map<int, boost::shared_ptr<Cube> >::const_iterator crvbars;
  for(crvbars=_crvscintillatorbars.begin(); crvbars!=_crvscintillatorbars.end(); crvbars++)
  {
    crvbars->second->makeGeometryVisible(visible);
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

  std::map<int,boost::shared_ptr<Cube> >::const_iterator crvbar;
  for(crvbar=_crvscintillatorbars.begin(); crvbar!=_crvscintillatorbars.end(); crvbar++)
  {
    crvbar->second->toForeground();
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
      if(color<=0 || std::isnan(color)) color=0;
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
      if(color<=0 || std::isnan(color)) color=0;
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
      if(color<=0 || std::isnan(color)) color=0;
      color+=2000;
      (*driftradius)->setColor(color);
    }
    else (*driftradius)->setColor(whitebackground?1:0);
  }
  std::vector<boost::shared_ptr<Cube> >::const_iterator crvhit;
  for(crvhit=_crvhits.begin(); crvhit!=_crvhits.end(); crvhit++)
  {
    double time=(*crvhit)->getStartTime();
    if(hitcolors)
    {
      int color=TMath::FloorNint(20.0*(time-mint)/(maxt-mint));
      if(color>=20) color=19;
      if(color<=0 || std::isnan(color)) color=0;
      color+=2000;
      (*crvhit)->setColor(color);
    }
    else (*crvhit)->setColor(whitebackground?1:0);
  }
}

void DataInterface::useTrackColors(boost::shared_ptr<ContentSelector> const &contentSelector, bool trackcolors, bool whitebackground)
{
  std::vector<ContentSelector::trackInfoStruct> selectedTracks=contentSelector->getSelectedTrackNames();
  TrackColorSelector colorSelector(&selectedTracks, whitebackground);
  std::vector<boost::shared_ptr<Track> >::const_iterator track;
  for(track=_tracks.begin(); track!=_tracks.end(); track++)
  {
    if(trackcolors)
    {
      int color=colorSelector.getColor(*track);
      (*track)->setColor(color);
    }
    else (*track)->setColor(whitebackground?kBlack:kWhite);
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

DataInterface::spaceminmax DataInterface::getSpaceBoundary(bool useTracker, bool useTarget, bool useCalorimeter, bool useTracks, bool useCRV)
{
  spaceminmax m;
  resetBoundaryP(m);
  if(useTracker)
  {
    findBoundaryP(m, _trackerMinmax.minx, _trackerMinmax.miny, _trackerMinmax.minz);
    findBoundaryP(m, _trackerMinmax.maxx, _trackerMinmax.maxy, _trackerMinmax.maxz);
  }
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
  if(useCRV)
  {
    findBoundaryP(m, _crvMinmax.minx, _crvMinmax.miny, _crvMinmax.minz);
    findBoundaryP(m, _crvMinmax.maxx, _crvMinmax.maxy, _crvMinmax.maxz);
  }

  if(std::isnan(m.minx)) m.minx=-1000;
  if(std::isnan(m.miny)) m.miny=-1000;
  if(std::isnan(m.minz)) m.minz=-1000;
  if(std::isnan(m.maxx)) m.maxx=1000;
  if(std::isnan(m.maxy)) m.maxy=1000;
  if(std::isnan(m.maxz)) m.maxz=1000;
  return m;
}

void DataInterface::findBoundaryT(timeminmax &m, double t)
{
  if(std::isnan(m.mint) || t<m.mint) m.mint=t;
  if(std::isnan(m.maxt) || t>m.maxt) m.maxt=t;
}

void DataInterface::findBoundaryP(spaceminmax &m, double x, double y, double z)
{
  if(std::isnan(m.minx) || x<m.minx) m.minx=x;
  if(std::isnan(m.miny) || y<m.miny) m.miny=y;
  if(std::isnan(m.minz) || z<m.minz) m.minz=z;
  if(std::isnan(m.maxx) || x>m.maxx) m.maxx=x;
  if(std::isnan(m.maxy) || y>m.maxy) m.maxy=y;
  if(std::isnan(m.maxz) || z>m.maxz) m.maxz=z;
}

using LHPT = mu2e::KalSeed::LHPT;
using CHPT = mu2e::KalSeed::CHPT;
using KLPT = mu2e::KalSeed::KLPT;

template<class KTRAJ> void DataInterface::fillKalSeedTrajectory(std::unique_ptr<KTRAJ> &trajectory,
                                                                int particleid, int trackclass, int trackclassindex, double p1,
                                                                boost::shared_ptr<ComponentInfo> info)
{
  double t1=trajectory->range().begin();
  double t2=trajectory->range().end();

  const auto &pos1=trajectory->position3(t1);
  const auto &pos2=trajectory->position3(t2);

  boost::shared_ptr<Track> track(new Track(pos1.x(),pos1.y(),pos1.z(),t1, pos2.x(),pos2.y(),pos2.z(),t2,
                                           particleid, trackclass, trackclassindex, p1,
                                           _geometrymanager, _topvolume, _mainframe, info, false));
  _components.push_back(track);
  _tracks.push_back(track);

  for(double t=t1; t<=t2; t+=_kalStepSize)
  {
    const auto &p = trajectory->position3(t);
    findBoundaryT(_tracksTimeMinmax, t);
    findBoundaryP(_tracksMinmax, p.x(), p.y(), p.z());
    track->addTrajectoryPoint(p.x(), p.y(), p.z(), t);
  }
}

void DataInterface::fillEvent(boost::shared_ptr<ContentSelector> const &contentSelector)
{
  auto const& ptable = mu2e::GlobalConstantsHandle<mu2e::ParticleDataList>();
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
      int sid = hit.strawId().asUint16();
      int trackid = hit.trackId().asInt();
      double time = hit.time();
      double energy = hit.eDep();
      std::map<int,boost::shared_ptr<Straw> >::iterator straw=_straws.find(sid);
      if(straw!=_straws.end() && !std::isnan(time))
      {
        double previousStartTime=straw->second->getStartTime();
        if(std::isnan(previousStartTime))
        {
          findBoundaryT(_hitsTimeMinmax, time);  //is it Ok to exclude all following hits from the time window?
          straw->second->setStartTime(time);
          straw->second->getComponentInfo()->setText(1,Form("hit time(s): %gns",time/CLHEP::ns));
          straw->second->getComponentInfo()->setText(2,Form("deposited energy(s): %geV",energy/CLHEP::eV));
          straw->second->getComponentInfo()->setText(3,Form("track id(s): %i",trackid));
          _hits.push_back(straw->second);
        }
        else
        {
          straw->second->getComponentInfo()->expandLine(1,Form("%gns",time/CLHEP::ns));
          straw->second->getComponentInfo()->expandLine(2,Form("%geV",energy/CLHEP::eV));
          straw->second->getComponentInfo()->expandLine(3,Form("%i",trackid));
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
      int sid  = hit.strawId().asUint16();
      double time = hit.time();
      double dt = hit.dt();
      double energy = hit.energyDep();
      std::map<int,boost::shared_ptr<Straw> >::iterator straw=_straws.find(sid);
      if(straw!=_straws.end() && !std::isnan(time))
      {
        double previousStartTime=straw->second->getStartTime();
        if(std::isnan(previousStartTime))
        {
          findBoundaryT(_hitsTimeMinmax, time);  //is it Ok to exclude all following hits from the time window?
          straw->second->setStartTime(time);
          straw->second->setHitNumber(hitnumber);
          straw->second->getComponentInfo()->setText(1,Form("hit time(s): %gns",time/CLHEP::ns));
          straw->second->getComponentInfo()->setText(2,Form("deposited energy(s): %geV",energy/CLHEP::eV));
          straw->second->getComponentInfo()->setText(3,Form("hit time interval(s): %gns",dt/CLHEP::ns));
          _hits.push_back(straw->second);
        }
        else
        {
          straw->second->getComponentInfo()->expandLine(1,Form("%gns",time/CLHEP::ns));
          straw->second->getComponentInfo()->expandLine(2,Form("%geV",energy/CLHEP::eV));
          straw->second->getComponentInfo()->expandLine(3,Form("%gns",dt/CLHEP::ns));
        }
      }
    }
  }
  /*
  const mu2e::KalRepCollection *kalRepHits=contentSelector->getSelectedHitCollection<mu2e::KalRepCollection>();
  if(kalRepHits!=nullptr)
  {
    boost::shared_ptr<TGraphErrors> residualGraph(new TGraphErrors());
    residualGraph->SetTitle("Residual Graph");

    for(unsigned int i=0; i<kalRepHits->size(); i++)
    {
      const KalRep &particle = kalRepHits->at(i);
      TrkHitVector const& hots = particle.hitVector();
      if(hots.size() > 0)
      {
        _numberHits+=hots.size();
        for(auto iter=hots.begin(); iter!=hots.end(); iter++)
        {
          const TrkHit *hitOnTrack = *iter;
          const mu2e::TrkStrawHit* strawHit = dynamic_cast<const mu2e::TrkStrawHit*>(hitOnTrack);
          if(strawHit)
          {
            int    sid=strawHit->straw().id().asUint16();
            double time = strawHit->time();
            double hitT0 = strawHit->hitT0()._t0; //this is the time the hit "arrived at the straw"
                                              //don't know what the other times are
            double strawtime = strawHit->comboHit().time();
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

            std::map<int,boost::shared_ptr<Straw> >::iterator straw=_straws.find(sid);
            if(straw!=_straws.end() && !std::isnan(time))
            {
              double previousStartTime=straw->second->getStartTime();
              if(std::isnan(previousStartTime))
              {
                findBoundaryT(_hitsTimeMinmax, hitT0);  //is it Ok to exclude all following hits from the time window?
                straw->second->setStartTime(hitT0);
                straw->second->getComponentInfo()->setText(1,Form("hitT0(s): %gns",hitT0/CLHEP::ns));
                straw->second->getComponentInfo()->setText(2,Form("hit time(s): %gns",time/CLHEP::ns));
                straw->second->getComponentInfo()->setText(3,Form("strawhit time(s): %gns",strawtime/CLHEP::ns));
                residualGraph->GetXaxis()->SetTitle("z [mm]");
                residualGraph->GetYaxis()->SetTitle("Residual [??]");
                straw->second->getComponentInfo()->getHistVector().push_back(boost::dynamic_pointer_cast<TObject>(residualGraph));
                _hits.push_back(straw->second);
              }
              else
              {
                straw->second->getComponentInfo()->expandLine(1,Form("%gns",hitT0/CLHEP::ns));
                straw->second->getComponentInfo()->expandLine(2,Form("%gns",time/CLHEP::ns));
                straw->second->getComponentInfo()->expandLine(3,Form("%gns",strawtime/CLHEP::ns));
              }

              const boost::shared_ptr<std::string> strawname=straw->second->getComponentInfo()->getName();
              boost::shared_ptr<ComponentInfo> info(new ComponentInfo());
              info->setName(Form("Drift Radius for %s",strawname->c_str()));
              info->setText(0,strawname->c_str());
              info->setText(1,Form("Drift Radius %gcm",driftRadius/CLHEP::cm));
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

  */
  // KalSeed hits
  const mu2e::KalSeedCollection *kalSeedsWithHits=contentSelector->getSelectedHitCollection<mu2e::KalSeedCollection>();
  if(kalSeedsWithHits!=NULL)
  {
    for(size_t i=0; i<kalSeedsWithHits->size(); i++)
    {
      const mu2e::KalSeed &kalseed = kalSeedsWithHits->at(i);
      const std::vector<mu2e::TrkStrawHitSeed> &hits = kalseed.hits();
      for(size_t j=0; j<hits.size(); j++)
      {
        const mu2e::TrkStrawHitSeed &hit = hits.at(j);
        int    sid = hit.strawId().asUint16();
        double time = hit.time();

        std::map<int,boost::shared_ptr<Straw> >::iterator straw=_straws.find(sid);
        if(straw!=_straws.end() && !std::isnan(time))
        {
          double previousStartTime=straw->second->getStartTime();
          if(std::isnan(previousStartTime))
          {
            findBoundaryT(_hitsTimeMinmax, time);  //is it Ok to exclude all following hits from the time window?
            straw->second->setStartTime(time);
            straw->second->getComponentInfo()->setText(1,Form("hit time(s): %gns",time/CLHEP::ns));
            _hits.push_back(straw->second);
          }
          else
          {
            straw->second->getComponentInfo()->expandLine(1,Form("%gns",time/CLHEP::ns));
          }
        }
      }
    }
  }

  // StepPoints at calorimeter
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
      if(crystal!=_crystals.end() && !std::isnan(time))
      {
        double previousStartTime=crystal->second->getStartTime();
        if(std::isnan(previousStartTime))
        {
          findBoundaryT(_hitsTimeMinmax, time);  //is it Ok to exclude all following hits from the time window?
          crystal->second->setStartTime(time);
          crystal->second->getComponentInfo()->setText(2,Form("hit time(s): %gns",time/CLHEP::ns));
          crystal->second->getComponentInfo()->setText(3,Form("deposited energy(s): %geV",energy/CLHEP::eV));
          crystal->second->getComponentInfo()->setText(4,Form("track ID(s): %i",trackid));
          _crystalhits.push_back(crystal->second);
        }
        else
        {
          crystal->second->getComponentInfo()->expandLine(2,Form("%gns",time/CLHEP::ns));
          crystal->second->getComponentInfo()->expandLine(3,Form("%geV",energy/CLHEP::eV));
          crystal->second->getComponentInfo()->expandLine(4,Form("%i",trackid));
        }
      }
    }
  }

  const mu2e::CaloHitCollection *calohits=contentSelector->getSelectedCaloHitCollection<mu2e::CaloHitCollection>();
  if(calohits!=nullptr)
  {
    _numberCrystalHits=calohits->size();  //this is not accurate since the return value gives the RO hits
    std::vector<mu2e::CaloHit>::const_iterator iter;
    for(iter=calohits->begin(); iter!=calohits->end(); iter++)
    {
      const mu2e::CaloHit& calohit = *iter;
      int crystalid = calohit.crystalID();

      double time = calohit.time();
      double energy = calohit.energyDep();
      std::map<int,boost::shared_ptr<VirtualShape> >::iterator crystal=_crystals.find(crystalid);
      if(crystal!=_crystals.end() && !std::isnan(time))
      {
        double previousStartTime=crystal->second->getStartTime();
        if(std::isnan(previousStartTime))
        {
          findBoundaryT(_hitsTimeMinmax, time);  //is it Ok to exclude all following hits from the time window?
          crystal->second->setStartTime(time);
          crystal->second->getComponentInfo()->setText(2,Form("hit time(s): %gns",time/CLHEP::ns));
          crystal->second->getComponentInfo()->setText(3,Form("deposited energy(s): %geV",energy/CLHEP::eV));
          _crystalhits.push_back(crystal->second);
        }
        else
        {
          crystal->second->getComponentInfo()->expandLine(2,Form("%gns",time/CLHEP::ns));
          crystal->second->getComponentInfo()->expandLine(3,Form("%geV",energy/CLHEP::eV));
        }
      }
    }
  }

//CRV waveforms
  for(std::map<int,boost::shared_ptr<Cube> >::iterator crvbar=_crvscintillatorbars.begin(); crvbar!=_crvscintillatorbars.end(); crvbar++)
  {
    crvbar->second->getComponentInfo()->getHistVector().clear();
  }

  double TDC0time = contentSelector->getTDC0time();
  const std::vector<art::Handle<mu2e::CrvDigiCollection> > &crvDigisVector = contentSelector->getSelectedCrvDigiCollection();
  for(size_t i=0; i<crvDigisVector.size(); i++)
  {
    const art::Handle<mu2e::CrvDigiCollection> &crvDigis = crvDigisVector[i];
    std::string moduleLabel = crvDigis.provenance()->moduleLabel();
    for(size_t j=0; j<crvDigis->size(); j++)
    {
      mu2e::CrvDigi const& digi(crvDigis->at(j));
      int index = digi.GetScintillatorBarIndex().asInt();
      int sipm  = digi.GetSiPMNumber();
      std::string multigraphName = Form("Waveform (%s) SiPM %i",moduleLabel.c_str(),sipm);
      std::map<int,boost::shared_ptr<Cube> >::iterator crvbar=_crvscintillatorbars.find(index);
      if(crvbar!=_crvscintillatorbars.end())
      {
        //each digi collection and each SiPM gets its own multigraph
        bool newMultigraph=true;
        int multigraphIndex=0;
        std::vector<boost::shared_ptr<TObject> > &v=crvbar->second->getComponentInfo()->getHistVector();
        for(size_t k=0; k<v.size(); k++)
        {
          //check whether multigraph exists already for this module label and SiPM
          if(multigraphName.compare(v[k]->GetName())==0) {newMultigraph=false; multigraphIndex=k; break;}
        }
        if(newMultigraph)
        {
          boost::shared_ptr<TMultiGraph> waveform(new TMultiGraph(multigraphName.c_str(),multigraphName.c_str()));
          v.push_back(boost::dynamic_pointer_cast<TObject>(waveform));
          multigraphIndex=v.size()-1;
        }

        TGraph *graph = new TGraph(digi.GetADCs().size());
        graph->SetMarkerStyle(20);
        graph->SetMarkerSize(2);
        for(size_t k=0; k<digi.GetADCs().size(); k++)
        {
          graph->SetPoint(k,TDC0time+(digi.GetStartTDC()+k)*mu2e::CRVDigitizationPeriod,digi.GetADCs()[k]);
        }
        boost::dynamic_pointer_cast<TMultiGraph>(v[multigraphIndex])->Add(graph,"p");
      }
    }
  }

//CRV reco pulses
  const mu2e::CrvRecoPulseCollection *crvRecoPulses=contentSelector->getSelectedCrvHitCollection<mu2e::CrvRecoPulseCollection>();
  if(crvRecoPulses!=nullptr)
  {
    for(size_t i=0; i<crvRecoPulses->size(); i++)
    {
      const mu2e::CrvRecoPulse &recoPulse = crvRecoPulses->at(i);
      int    index = recoPulse.GetScintillatorBarIndex().asInt();
      int    sipm  = recoPulse.GetSiPMNumber();
      double time  = recoPulse.GetPulseTime();
      int    PEs   = recoPulse.GetPEs();

      size_t channel  = index*4 + sipm;
      double recoPulsePedestal = _calib->pedestal(channel);

      std::map<int,boost::shared_ptr<Cube> >::iterator crvbar=_crvscintillatorbars.find(index);
      if(crvbar!=_crvscintillatorbars.end() && !std::isnan(time))
      {
        double previousStartTime=crvbar->second->getStartTime();
        if(std::isnan(previousStartTime)) _crvhits.push_back(crvbar->second);  //first reco hit of this counter

        if(std::isnan(previousStartTime) || time<previousStartTime) crvbar->second->setStartTime(time);

        if(std::isnan(previousStartTime))
        {
          findBoundaryT(_hitsTimeMinmax, time);  //is it Ok to exclude all following hits from the time window?
          crvbar->second->getComponentInfo()->setText(3,"Reco pulse SiPM0 PEs/time: ");
          crvbar->second->getComponentInfo()->setText(4,"Reco pulse SiPM1 PEs/time: ");
          crvbar->second->getComponentInfo()->setText(5,"Reco pulse SiPM2 PEs/time: ");
          crvbar->second->getComponentInfo()->setText(6,"Reco pulse SiPM3 PEs/time: ");
        }
        crvbar->second->getComponentInfo()->expandLine(sipm+3,Form("%iPEs/%.1fns",PEs,time/CLHEP::ns));

        //each digi collection and each SiPM gets its own multigraph
        std::vector<boost::shared_ptr<TObject> > &v=crvbar->second->getComponentInfo()->getHistVector();
        for(size_t k=0; k<v.size(); k++)
        {
          //check whether this multigraph is for the current SiPM
          const char *multigraphName = v[k]->GetName();
          int nameLength = strlen(multigraphName);
          if(nameLength<1) continue;
          if(atoi(multigraphName+nameLength-1)==sipm)
          {
            double fitParam0 = recoPulse.GetPulseHeight()*TMath::E();
            double fitParam1 = recoPulse.GetPulseTime();
            double fitParam2 = recoPulse.GetPulseBeta();

            TList *functionList = boost::dynamic_pointer_cast<TMultiGraph>(v[k])->GetListOfFunctions();
            if(functionList->GetSize()==0) functionList->Add(new TF1("pedestal",Form("%f",recoPulsePedestal))); //TODO: Use compiled function
            TF1 *f = new TF1(Form("peakfitter%i",functionList->GetSize()),
                             Form("%f*(TMath::Exp(-(x-%f)/%f-TMath::Exp(-(x-%f)/%f)))",
                             fitParam0,fitParam1,fitParam2,fitParam1,fitParam2));  //TODO: Use compiled function
            f->SetLineWidth(2);
            f->SetLineColor(2);
            if(recoPulse.GetRecoPulseFlags().test(mu2e::CrvRecoPulseFlagEnums::failedFit)) f->SetLineStyle(2);
            functionList->Add(f);
          }
        }
      }
    }
  }

  const mu2e::PhysicalVolumeInfoMultiCollection *physicalVolumesMulti=contentSelector->getPhysicalVolumeInfoMultiCollection();

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
      art::Ptr<mu2e::SimParticle> particlePtr(trackInfos[i].productId, &particle, iter->first.asUint());
      int id = particle.id().asInt();   //is identical with cet::map_vector_key& particleKey = iter->first;
      int parentid = particle.parentId().asInt();
      int particleid=particle.pdgId();
      int trackclass=trackInfos[i].classID;
      int trackclassindex=trackInfos[i].index;
      std::string particlecollection=trackInfos[i].entryText;
      std::string particlename=ptable->particle(particle.pdgId()).name();
      std::string startVolumeName="unknown volume";
      std::string endVolumeName="unknown volume";
      if(physicalVolumesMulti!=nullptr)
      {
        mu2e::PhysicalVolumeMultiHelper volumeMultiHelper(physicalVolumesMulti);
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

      boost::shared_ptr<ComponentInfo> info(new ComponentInfo());
      info->setName(Form("Track %i  %s  (%s)",id,particlename.c_str(),particlecollection.c_str()));
      if(parentid!=0) info->setText(0,Form("Track %i  %s  Parent %i  (%s)",id,particlename.c_str(),parentid,particlecollection.c_str()));
      else info->setText(0,Form("Track %i  %s  generated track  (%s)",id,particlename.c_str(),particlecollection.c_str()));
      info->setText(1,Form("Start Energy %gMeV  End Energy %gMeV",e1/CLHEP::MeV,e2/CLHEP::MeV));
      info->setText(2,Form("Created by %s in %s",particle.creationCode().name().c_str(),startVolumeName.c_str()));
      info->setText(3,Form("Destroyed by %s in %s",particle.stoppingCode().name().c_str(),endVolumeName.c_str()));
      info->setText(4,"Daughter IDs:");
      std::vector<art::Ptr<mu2e::SimParticle> >::const_iterator daughter;
      for(daughter=particle.daughters().begin();
          daughter!=particle.daughters().end();
          daughter++)
      {
        info->expandLine(4,Form("%lu",(*daughter)->id().asInt()));
      }
      boost::shared_ptr<Track> shape(new Track(x1,y1,z1,t1, x2,y2,z2,t2,
                                               particleid, trackclass, trackclassindex, e1,
                                               _geometrymanager, _topvolume, _mainframe, info, false));
      findTrajectory(contentSelector,shape,particle.id(), trackInfos[i]);
      _components.push_back(shape);
      _tracks.push_back(shape);
    }
  }
  /*
  trackInfos.clear();
  std::vector<const mu2e::KalRepCollection*> kalRepCollectionVector=contentSelector->getSelectedTrackCollection<mu2e::KalRepCollection>(trackInfos);
  for(unsigned int i=0; i<kalRepCollectionVector.size(); i++)
  {
    const mu2e::KalRepCollection *kalReps=kalRepCollectionVector[i];
    for(unsigned int j=0; j<kalReps->size(); j++)
    {
      KalRep const* kalrep = kalReps->get(j);
        int trackclass=trackInfos[i].classID;
        int trackclassindex=trackInfos[i].index;
        std::string trackcollection=trackInfos[i].entryText;
        int particleid=kalrep->particleType().particleType();
        std::string particlename=ptable->particle(particleid).name();
        std::string c=Form("Kalman Track %i  %s  (%s)",j,particlename.c_str(),trackcollection.c_str());
        boost::shared_ptr<ComponentInfo> info(new ComponentInfo());
        info->setName(c.c_str());
        info->setText(0,c.c_str());

        double hitcount=0;
        double offset=0;
        TrkHitVector const& hots=kalrep->hitVector();
        if(hots.size()>0)
        {
          boost::shared_ptr<TGraphErrors> residualGraph(new TGraphErrors());
          residualGraph->SetTitle("Residual Graph");
          for(auto iter=hots.begin(); iter!=hots.end(); iter++)
          {
            const TrkHit *hitOnTrack = *iter;
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
        double t0=kalrep->t0().t0();
        double firsthitfltlen = kalrep->lowFitRange();
        double lasthitfltlen = kalrep->hiFitRange();
        double entlen = min(firsthitfltlen,lasthitfltlen);
        double loclen(0.0);
        const TrkSimpTraj* ltraj = kalrep->localTrajectory(entlen,loclen);
        const CLHEP::HepVector &params=ltraj->parameters()->parameter();
        double d0 = params[0];
        double om = params[2];
        double rmax = d0+2.0/om;

        info->setText(1,Form("Charge %i",charge));
        info->setText(2,Form("Start Momentum %gMeV/c  End Momentum %gMeV/c",p1/CLHEP::MeV,p2/CLHEP::MeV));
        info->setText(3,Form("t0 %gns  d0 %gmm  rmax %gmm",t0/CLHEP::ns,d0/CLHEP::mm,rmax/CLHEP::mm));
    }
  }
  */
  // KalSeed tracks
  trackInfos.clear();
  std::vector<const mu2e::KalSeedCollection*> kalSeedCollectionVector=contentSelector->getSelectedTrackCollection<mu2e::KalSeedCollection>(trackInfos);
  for(size_t i=0; i<kalSeedCollectionVector.size(); i++)
  {
    const mu2e::KalSeedCollection *kalSeeds=kalSeedCollectionVector[i];
    for(size_t j=0; j<kalSeeds->size(); j++)
    {
      const mu2e::KalSeed &kalseed = kalSeeds->at(j);
      int trackclass=trackInfos[i].classID;
      int trackclassindex=trackInfos[i].index;
      std::string trackcollection=trackInfos[i].entryText;
      int particleid=kalseed.particle();
      std::string particlename=ptable->particle(particleid).name();

      const std::vector<mu2e::KalSegment> &segments = kalseed.segments();
      size_t nSegments=segments.size();
      if(nSegments==0) continue;
      const mu2e::KalSegment &segmentFirst = kalseed.segments().front();
      const mu2e::KalSegment &segmentLast = kalseed.segments().back();
      double fltLMin=segmentFirst.fmin();
      double fltLMax=segmentLast.fmax();
      XYZVectorF momvec1, momvec2;
      segmentFirst.mom(fltLMin, momvec1);
      segmentLast.mom(fltLMax, momvec2);
      double p1=mu2e::GenVector::Hep3Vec(momvec1).mag();
      double p2=mu2e::GenVector::Hep3Vec(momvec2).mag();

      boost::shared_ptr<ComponentInfo> info(new ComponentInfo());
      std::string c=Form("KalSeed Track %lu  %s  (%s)",j,particlename.c_str(),trackcollection.c_str());
      info->setName(c.c_str());
      info->setText(0,c.c_str());
      info->setText(1,Form("Start Momentum %gMeV/c  End Momentum %gMeV/c",p1/CLHEP::MeV,p2/CLHEP::MeV));

      if(kalseed.loopHelixFit())
      {
        info->setText(2,"Loop Helix Fit");
        auto trajectory=kalseed.loopHelixFitTrajectory();
        fillKalSeedTrajectory<LHPT>(trajectory, particleid, trackclass, trackclassindex, p1, info);
      }
      else if(kalseed.centralHelixFit())
      {
        info->setText(2,"Central Helix Fit");
        auto trajectory=kalseed.centralHelixFitTrajectory();
        fillKalSeedTrajectory<CHPT>(trajectory, particleid, trackclass, trackclassindex, p1, info);
      }
      else if(kalseed.kinematicLineFit())
      {
        info->setText(2,"Kinematic Line Fit");
        auto trajectory=kalseed.kinematicLineFitTrajectory();
        fillKalSeedTrajectory<KLPT>(trajectory, particleid, trackclass, trackclassindex, p1, info);
      }
      else
      {
//old KalSeed
        info->setText(2,"Old KalSeed using Segments");

        double t0   = kalseed.t0().t0();
        double flt0 = kalseed.flt0();
        double mass = ptable->particle(kalseed.particle()).mass();
        double v  = mu2e::TrkUtilities::beta(mass,(p1+p2)/2.0)*CLHEP::c_light;

        for(size_t k=0; k<nSegments; k++)
        {
          const mu2e::KalSegment &segment = segments.at(k);

          fltLMin=segment.fmin();
          fltLMax=segment.fmax();

/*
//interpolation between segments doesn't seem to work anymore
          if(k>0)
          {
            double fltLMaxPrev=segments.at(k-1).fmax();
            fltLMin=(fltLMin+fltLMaxPrev)/2.0;
          }
          if(k+1<nSegments)
          {
            double fltLMinNext=segments.at(k+1).fmin();
            fltLMax=(fltLMax+fltLMinNext)/2.0;
          }
*/

          XYZVectorF pos1, pos2;
          segment.helix().position(fltLMin,pos1);
          segment.helix().position(fltLMax,pos2);
          double x1=mu2e::GenVector::Hep3Vec(pos1).x();
          double y1=mu2e::GenVector::Hep3Vec(pos1).y();
          double z1=mu2e::GenVector::Hep3Vec(pos1).z();
          double x2=mu2e::GenVector::Hep3Vec(pos2).x();
          double y2=mu2e::GenVector::Hep3Vec(pos2).y();
          double z2=mu2e::GenVector::Hep3Vec(pos2).z();
          double t1=t0+(fltLMin-flt0)/v;
          double t2=t0+(fltLMax-flt0)/v;
          boost::shared_ptr<Track> track(new Track(x1,y1,z1,t1, x2,y2,z2,t2,
                                                   particleid, trackclass, trackclassindex, p1,
                                                   _geometrymanager, _topvolume, _mainframe, info, false));
          _components.push_back(track);
          _tracks.push_back(track);

          for(double fltL=fltLMin; fltL<=fltLMax; fltL+=1.0)
          {
            double t=t0+(fltL-flt0)/v;
            XYZVectorF pos;
            segment.helix().position(fltL,pos);
            CLHEP::Hep3Vector p = mu2e::GenVector::Hep3Vec(pos);
            findBoundaryT(_tracksTimeMinmax, t);
            findBoundaryP(_tracksMinmax, p.x(), p.y(), p.z());
            track->addTrajectoryPoint(p.x(), p.y(), p.z(), t);
          }
        } //segments
      } //old KalSeed
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

      std::string particlename=ptable->particle(particleid).name();
      boost::shared_ptr<ComponentInfo> info(new ComponentInfo());
      std::string c=Form("TrkExt Trajectory %i  %s  (%s)", trkExtTraj.id(), particlename.c_str(),trackcollection.c_str());
      info->setName(c.c_str());
      info->setText(0,c.c_str());

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
}

void DataInterface::findTrajectory(boost::shared_ptr<ContentSelector> const &contentSelector,
                                   boost::shared_ptr<Track> const &track, const cet::map_vector_key &id,
                                   const ContentSelector::trackInfoStruct &trackInfo)
{
  const mu2e::MCTrajectoryCollection *mcTrajectories=contentSelector->getMCTrajectoryCollection(trackInfo);
  if(mcTrajectories!=nullptr)
  {
    std::map<art::Ptr<mu2e::SimParticle>,mu2e::MCTrajectory>::const_iterator traj_iter;
    for(traj_iter=mcTrajectories->begin(); traj_iter!=mcTrajectories->end(); traj_iter++)
    {
      if(traj_iter->first.key()==id.asUint())
//      if(traj_iter->first->id()==id)
//      if(traj_iter->second.sim()->id()==id)
      {
        const auto& points = traj_iter->second.points();
        for(auto point_iter=points.begin(); point_iter!=points.end(); ++point_iter)
        {
          track->addTrajectoryPoint(point_iter->x()-_detSysOrigin.x(),
                                    point_iter->y()-_detSysOrigin.y(),
                                    point_iter->z()-_detSysOrigin.z(),
                                    point_iter->t());
        }
      }
    }
  }
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
    for(int i=1; i<7; i++) (*hit)->getComponentInfo()->removeLine(i);  //keep first line
    (*hit)->getComponentInfo()->getHistVector().clear();
    (*hit)->setHitNumber(-1);
    (*hit)->setStartTime(NAN);
    (*hit)->start();
  }
  std::vector<boost::shared_ptr<VirtualShape> >::const_iterator crystalhit;
  for(crystalhit=_crystalhits.begin(); crystalhit!=_crystalhits.end(); crystalhit++)
  {
    for(int i=1; i<7; i++) (*crystalhit)->getComponentInfo()->removeLine(i); //keep first line
    (*crystalhit)->getComponentInfo()->getHistVector().clear();
    (*crystalhit)->setStartTime(NAN);
    (*crystalhit)->start();
  }
  std::vector<boost::shared_ptr<Cube> >::const_iterator crvhit;
  for(crvhit=_crvhits.begin(); crvhit!=_crvhits.end(); crvhit++)
  {
    for(int i=1; i<7; i++) (*crvhit)->getComponentInfo()->removeLine(i); //keep first line
    (*crvhit)->getComponentInfo()->getHistVector().clear();
    (*crvhit)->setStartTime(NAN);
    (*crvhit)->start();
  }

  _hits.clear();
  _crystalhits.clear();
  _crvhits.clear();
  _tracks.clear();  //will call the d'tors of all tracks, since they aren't used anywhere anymore
  _driftradii.clear(); //will call the d'tors of all driftradii, since they aren't used anywhere anymore

  _mainframe->getHistDrawVector().clear();
}

void DataInterface::removeAllComponents()
{
  _components.clear();
  _straws.clear();
  _crystals.clear();
  _crvscintillatorbars.clear();
  _hits.clear();
  _crystalhits.clear();
  _crvhits.clear();
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
