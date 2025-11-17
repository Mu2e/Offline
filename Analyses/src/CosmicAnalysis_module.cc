#include <iostream>
#include <memory>
#include <string>

#include <TDirectory.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>

#include "CLHEP/Units/GlobalSystemOfUnits.h"
#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "Offline/DataProducts/inc/VirtualDetectorId.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/GeometryService/inc/VirtualDetector.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/PhysicalVolumeInfoMultiCollection.hh"
#include "Offline/MCDataProducts/inc/MCTrajectoryCollection.hh"
#include "Offline/MCDataProducts/inc/PtrStepPointMCVector.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/Mu2eUtilities/inc/PhysicalVolumeMultiHelper.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/CrvCoincidence.hh"
#include "Offline/RecoDataProducts/inc/CrvRecoPulse.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/ParameterSet.h"

using namespace CLHEP;
#include "Offline/RecoDataProducts/inc/KalRepCollection.hh"

typedef struct
{
  double cosmic_pos[3];
  double cosmic_E, cosmic_p;
  double cosmic_ph, cosmic_th, cosmic_costh;
  char   cosmic_particle[100];

  int    reco_n; //number of reconstructed tracks
  double reco_t0, reco_p0;

  char   simreco_particle[100];
  char   simreco_production_process[100], simreco_production_volume[100];
  double simreco_startp, simreco_endp, simreco_startp_z;
  double simreco_pos[3];

  double xplane1[3], xplane2[3], xplane3[3];
  double yplane1[3];
  double zplane1[3], zplane2[3], zplane3[3];

  double xplane1Dir[3], xplane2Dir[3], xplane3Dir[3];
  double yplane1Dir[3];
  double zplane1Dir[3], zplane2Dir[3], zplane3Dir[3];

  double firstCoincidenceHitTime;
  int    firstCoincidenceHitSectorType;
  double firstCoincidenceHitPos[3];

  bool   CRVveto_allSectors;
  bool   CRVveto[8];
  double CRVvetoTime[8];
  double CRVvetoPos[8][3];

  bool   CRVhit_allSectors;
  bool   CRVhit[8];
  double CRVhitTime[8];
  double CRVhitPos[8][3];
  double CRVhitDir[8][3];

  int    run_number, subrun_number, event_number;

  char filename[200];

  void clear()
  {
    reco_n=0;
    reco_t0=NAN;
    reco_p0=NAN;
    filename[0]='\0';
    simreco_particle[0]='\0';
    simreco_production_process[0]='\0';
    simreco_production_volume[0]='\0';
    simreco_startp=NAN;
    simreco_endp=NAN;
    simreco_startp_z=NAN;

    for(int i=0; i<3; i++)
    {
      simreco_pos[i]=NAN;
      xplane1[i]=NAN;
      xplane2[i]=NAN;
      xplane3[i]=NAN;
      yplane1[i]=NAN;
      zplane1[i]=NAN;
      zplane2[i]=NAN;
      zplane3[i]=NAN;
      xplane1Dir[i]=NAN;
      xplane2Dir[i]=NAN;
      xplane3Dir[i]=NAN;
      yplane1Dir[i]=NAN;
      zplane1Dir[i]=NAN;
      zplane2Dir[i]=NAN;
      zplane3Dir[i]=NAN;
      firstCoincidenceHitPos[i]=NAN;
    }

    firstCoincidenceHitTime=NAN;

    CRVveto_allSectors=false;
    CRVhit_allSectors=false;
    for(int i=0; i<8; i++)
    {
      CRVveto[i]=false;
      CRVvetoTime[i]=NAN;
      CRVhit[i]=false;
      CRVhitTime[i]=NAN;
      for(int j=0; j<3; j++)
      {
        CRVvetoPos[i][j]=NAN;
        CRVhitPos[i][j]=NAN;
        CRVhitDir[i][j]=NAN;
      }
    }
  };
} EventInfo;

namespace mu2e
{
  class CosmicAnalysis : public art::EDAnalyzer
  {
    public:
    explicit CosmicAnalysis(fhicl::ParameterSet const &pset);
    virtual ~CosmicAnalysis() { }
    virtual void beginJob();
    void analyze(const art::Event& e);
    void findCrossingDetails(const std::vector<MCTrajectoryPoint> &trajectoryPoints, int dim, double crossingPos,
                             double *crossingPoint, double *crossingDirection);
    void findCrossings(const art::Event& event, const cet::map_vector_key& particleKey);

    private:
    std::string _fitterModuleLabel;
    std::string _fitterModuleInstance;
    std::string _g4ModuleLabel;
    std::string _generatorModuleLabel;
    std::string _hitmakerModuleLabel;
    std::string _hitmakerModuleInstance;
    std::string _crvCoincidenceModuleLabel;
    std::string _crvRecoPulseModuleLabel;
    std::string _volumeModuleLabel;
    std::string _filenameSearchPattern;
    EventInfo   _eventinfo;
    TTree      *_tree;

    CLHEP::Hep3Vector _detSysOrigin;
    double            _zent;
  };

  CosmicAnalysis::CosmicAnalysis(fhicl::ParameterSet const &pset)
    :
    art::EDAnalyzer(pset),
    _fitterModuleLabel(pset.get<std::string>("fitterModuleLabel")),
    _fitterModuleInstance(pset.get<std::string>("fitterModuleInstance")),
    _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel","detectorFilter")),
    _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel")),
    _hitmakerModuleLabel(pset.get<std::string>("hitmakerModuleLabel")),
    _hitmakerModuleInstance(pset.get<std::string>("hitmakerModuleInstance")),
    _crvCoincidenceModuleLabel(pset.get<std::string>("crvCoincidenceModuleLabel")),
    _crvRecoPulseModuleLabel(pset.get<std::string>("crvRecoPulseModuleLabel")),
    _volumeModuleLabel(pset.get<std::string>("volumeModuleLabel")),
    _filenameSearchPattern(pset.get<std::string>("filenameSearchPattern"))
  {
  }

  void CosmicAnalysis::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory tfdir = tfs->mkdir("cosmicInfo");
    _tree = tfdir.make<TTree>("cosmicInfoTree","cosmicInfoTree");

    EventInfo &e=_eventinfo;
    _tree->Branch("cosmic_pos",e.cosmic_pos,"cosmic_pos[3]/D");
    _tree->Branch("cosmic_E",&e.cosmic_E,"cosmic_E/D");
    _tree->Branch("cosmic_p",&e.cosmic_p,"cosmic_p/D");
    _tree->Branch("cosmic_ph",&e.cosmic_ph,"cosmic_ph/D");
    _tree->Branch("cosmic_th",&e.cosmic_th,"cosmic_th/D");
    _tree->Branch("cosmic_costh",&e.cosmic_costh,"cosmic_costh/D");
    _tree->Branch("cosmic_particle",e.cosmic_particle,"cosmic_particle[100]/C");
    _tree->Branch("reco_n",&e.reco_n,"reco_n/I");
    _tree->Branch("reco_t0",&e.reco_t0,"reco_t0/D");
    _tree->Branch("reco_p0",&e.reco_p0,"reco_p0/D");
    _tree->Branch("simreco_particle",e.simreco_particle,"simreco_particle[100]/C");
    _tree->Branch("simreco_production_process",e.simreco_production_process,"simreco_production_process[100]/C");
    _tree->Branch("simreco_production_volume",e.simreco_production_volume,"simreco_production_volume[100]/C");
    _tree->Branch("simreco_startp",&e.simreco_startp,"simreco_startp/D");
    _tree->Branch("simreco_endp",&e.simreco_endp,"simreco_endp/D");
    _tree->Branch("simreco_startp_z",&e.simreco_startp_z,"simreco_startp_z/D");
    _tree->Branch("simreco_pos",e.simreco_pos,"simreco_pos[3]/D");
    _tree->Branch("xplane1",e.xplane1,"xplane1[3]/D");
    _tree->Branch("xplane2",e.xplane2,"xplane2[3]/D");
    _tree->Branch("xplane3",e.xplane3,"xplane3[3]/D");
    _tree->Branch("yplane1",e.yplane1,"yplane1[3]/D");
    _tree->Branch("zplane1",e.zplane1,"zplane1[3]/D");
    _tree->Branch("zplane2",e.zplane2,"zplane2[3]/D");
    _tree->Branch("zplane3",e.zplane3,"zplane3[3]/D");
    _tree->Branch("xplane1Dir",e.xplane1Dir,"xplane1Dir[3]/D");
    _tree->Branch("xplane2Dir",e.xplane2Dir,"xplane2Dir[3]/D");
    _tree->Branch("xplane3Dir",e.xplane3Dir,"xplane3Dir[3]/D");
    _tree->Branch("yplane1Dir",e.yplane1Dir,"yplane1Dir[3]/D");
    _tree->Branch("zplane1Dir",e.zplane1Dir,"zplane1Dir[3]/D");
    _tree->Branch("zplane2Dir",e.zplane2Dir,"zplane2Dir[3]/D");
    _tree->Branch("zplane3Dir",e.zplane3Dir,"zplane3Dir[3]/D");
    _tree->Branch("firstCoincidenceHitTime",&e.firstCoincidenceHitTime,"firstCoincidenceHitTime/D");
    _tree->Branch("firstCoincidenceHitSectorType",&e.firstCoincidenceHitSectorType,"firstCoincidenceHitSectorType/I");
    _tree->Branch("firstCoincidenceHitPos",e.firstCoincidenceHitPos,"firstCoincidenceHitPos[3]/D");
    _tree->Branch("CRVveto_allSectors",&e.CRVveto_allSectors,"CRVveto_allSectors/O");
    _tree->Branch("CRVveto",e.CRVveto,"CRVveto[8]/O");
    _tree->Branch("CRVvetoTime",e.CRVvetoTime,"CRVvetoTime[8]/D");
    _tree->Branch("CRVvetoPos",e.CRVvetoPos,"CRVvetoPos[8][3]/D");
    _tree->Branch("CRVhit_allSectors",&e.CRVhit_allSectors,"CRVhit_allSectors/O");
    _tree->Branch("CRVhit",e.CRVhit,"CRVhit[8]/O");
    _tree->Branch("CRVhitTime",e.CRVhitTime,"CRVhitTime[8]/D");
    _tree->Branch("CRVhitPos",e.CRVhitPos,"CRVhitPos[8][3]/D");
    _tree->Branch("CRVhitDir",e.CRVhitDir,"CRVhitDir[8][3]/D");
    _tree->Branch("run_number",&e.run_number,"run_number/I");
    _tree->Branch("subrun_number",&e.subrun_number,"subrun_number/I");
    _tree->Branch("event_number",&e.event_number,"event_number/I");
    _tree->Branch("filename",e.filename,"filename[200]/C");
  }

  void CosmicAnalysis::findCrossingDetails(const std::vector<MCTrajectoryPoint> &trajectoryPoints, int dim, double crossingPos,
                                           double *crossingPoint, double *crossingDirection)
  {
    if(!std::isnan(crossingPoint[0])) return;  //point already found
    for(unsigned int i=1; i<trajectoryPoints.size(); i++)
    {
      CLHEP::Hep3Vector point1=trajectoryPoints[i-1].pos()-_detSysOrigin;
      CLHEP::Hep3Vector point2=trajectoryPoints[i].pos()-_detSysOrigin;
      CLHEP::Hep3Vector diffVector=point2-point1;
      if(diffVector[dim]==0) continue;  //these two points are both on the same plane, try to find another pair
      if((point1[dim]>=crossingPos && point2[dim]<=crossingPos)
      || (point1[dim]<=crossingPos && point2[dim]>=crossingPos))
      {
        double ratio=(crossingPos - point1[dim])/diffVector[dim];
        CLHEP::Hep3Vector point=ratio*diffVector+point1;
        for(int i=0; i<3; i++)
        {
          crossingPoint[i]=point[i];
          crossingDirection[i]=diffVector.unit()[i];
        }
        return;  //don't look for a second point
      }
    }
  }

  void CosmicAnalysis::findCrossings(const art::Event& event, const cet::map_vector_key& particleKey)
  {
    double xCrossing1= -2525.4; //CRV-R
    double xCrossing2=  2525.4; //CRV-L
    double xCrossing3=  3800.0; //TS entrance
    double yCrossing1=  2589.7; //CRV-T / CRV-TS
    double zCrossing1=-12496.8; //CRV-U
    double zCrossing2= -7512.6; //open left side at the TS entrance
    double zCrossing3=  8571.0; //CRV-D

    std::string productInstanceName="";
    art::Handle<mu2e::MCTrajectoryCollection> _mcTrajectories;
    if(event.getByLabel(_g4ModuleLabel,productInstanceName,_mcTrajectories))
    {
      std::map<art::Ptr<mu2e::SimParticle>,mu2e::MCTrajectory>::const_iterator traj_iter;
      for(traj_iter=_mcTrajectories->begin(); traj_iter!=_mcTrajectories->end(); traj_iter++)
      {
        if(traj_iter->first->id()==particleKey)
        {
          const auto &trajectoryPoints = traj_iter->second.points();
          findCrossingDetails(trajectoryPoints, 0, xCrossing1, _eventinfo.xplane1, _eventinfo.xplane1Dir);
          findCrossingDetails(trajectoryPoints, 0, xCrossing2, _eventinfo.xplane2, _eventinfo.xplane2Dir);
          findCrossingDetails(trajectoryPoints, 0, xCrossing3, _eventinfo.xplane3, _eventinfo.xplane3Dir);
          findCrossingDetails(trajectoryPoints, 1, yCrossing1, _eventinfo.yplane1, _eventinfo.yplane1Dir);
          findCrossingDetails(trajectoryPoints, 2, zCrossing1, _eventinfo.zplane1, _eventinfo.zplane1Dir);
          findCrossingDetails(trajectoryPoints, 2, zCrossing2, _eventinfo.zplane2, _eventinfo.zplane2Dir);
          findCrossingDetails(trajectoryPoints, 2, zCrossing3, _eventinfo.zplane3, _eventinfo.zplane3Dir);
        }
      }
    }
  }

  inline bool operator < (const art::Ptr<mu2e::SimParticle> l, const art::Ptr<mu2e::SimParticle> r)
  {
    return l->id()<r->id();
  }

  void CosmicAnalysis::analyze(const art::Event& event)
  {
//FIXME
//Will be fixed later after the update of the CRV code is done
/*
    _eventinfo.clear();

    for(int i=0; i<gROOT->GetListOfFiles()->GetEntries(); i++)
    {
      std::string filename=((TFile*)gROOT->GetListOfFiles()->At(i))->GetName();
      if(filename.find(_filenameSearchPattern)!=std::string::npos) strncpy(_eventinfo.filename, filename.c_str(),199);
    }

    _eventinfo.run_number=event.id().run();
    _eventinfo.subrun_number=event.id().subRun();
    _eventinfo.event_number=event.id().event();

    GeomHandle<DetectorSystem> det;
    GeomHandle<VirtualDetector> vdg;
    _detSysOrigin = det->getOrigin();
    _zent = det->toDetector(vdg->getGlobal(VirtualDetectorId::TT_FrontPA)).z();

//generated cosmic ray muons
    art::Handle<GenParticleCollection> genParticles;
    if(event.getByLabel(_generatorModuleLabel, genParticles))
    {
      if(genParticles.product()!=NULL)
      {
        if(genParticles->size()>0)   //there shouldn't be more than 1 particle
        {
          const mu2e::GenParticle &particle=genParticles->at(0);
          _eventinfo.cosmic_pos[0]=particle.position().x()-_detSysOrigin.x();
          _eventinfo.cosmic_pos[1]=particle.position().y()-_detSysOrigin.y();
          _eventinfo.cosmic_pos[2]=particle.position().z()-_detSysOrigin.z();
          _eventinfo.cosmic_E=particle.momentum().e();

          double p=particle.momentum().vect().mag();
          double px=particle.momentum().px();
          double py=particle.momentum().py();
          double pz=particle.momentum().pz();
          double zenith=acos(py/p)*180.0/TMath::Pi();
          double azimuth=atan2(px,pz)*180.0/TMath::Pi();
          if(azimuth<0) azimuth+=360.0;

          _eventinfo.cosmic_p=p;
          _eventinfo.cosmic_ph=azimuth;
          _eventinfo.cosmic_th=zenith;
          _eventinfo.cosmic_costh=py/p;

          PDGCode::type pdgId = particle.pdgId();
          strncpy(_eventinfo.cosmic_particle, HepPID::particleName(pdgId).c_str(),99);
        }
      }
    }

//particles crossing the CRV planes
    std::string productInstanceName="";
    art::Handle<SimParticleCollection> simParticleCollection;
    if(event.getByLabel(_g4ModuleLabel,productInstanceName,simParticleCollection))
    {
      if(simParticleCollection.product()!=NULL)
      {
        cet::map_vector<mu2e::SimParticle>::const_iterator iter;
        for(iter=simParticleCollection->begin(); iter!=simParticleCollection->end(); iter++)
        {
          const mu2e::SimParticle& particle = iter->second;
          int pdgId=particle.pdgId();
          if(abs(pdgId)==11 || abs(pdgId)==13)
          {
            const cet::map_vector_key& particleKey = iter->first;
            findCrossings(event,particleKey); //if a second particle crosses the same plane, it will be ignored
          }
        }
      }
    }

//get the physical volume information for later
    art::Handle<mu2e::PhysicalVolumeInfoMultiCollection> physicalVolumesMulti;
    bool hasPhysicalVolumes=event.getSubRun().getByLabel(_volumeModuleLabel, physicalVolumesMulti);

//find the straw hits used by the reconstructed track
    std::map<art::Ptr<mu2e::SimParticle>, int> simParticles;   // < operator defined above

    art::Handle<KalRepCollection> kalReps;
    art::Handle<PtrStepPointMCVectorCollection> stepPointVectors;
    if(event.getByLabel(_fitterModuleLabel,_fitterModuleInstance,kalReps) && event.getByLabel(_hitmakerModuleLabel,_hitmakerModuleInstance,stepPointVectors))
    {
      if(kalReps.product()!=NULL && stepPointVectors.product()!=NULL)
      {
        if(kalReps->size()>0)
        {

          //if multiple reco tracks, selected reco track which has a momentum which is closest to 104.375 MeV/c
          size_t selectedTrack=0;
          double minMomentumDifference=NAN;
          for(size_t k=0; k<kalReps->size(); k++)
          {
            double momentumDifference=fabs(kalReps->at(k).momentum(0).mag()-104.375);
            if(momentumDifference<minMomentumDifference || std::isnan(minMomentumDifference)) selectedTrack=k;
          }

          const KalRep &particle = kalReps->at(selectedTrack);
          _eventinfo.reco_n=kalReps->size();
          _eventinfo.reco_t0=particle.t0().t0();

          //from TrkDiag/src/KalDiag.cc
          double firsthitfltlen = particle.lowFitRange();
          double lasthitfltlen = particle.hiFitRange();
          double entlen = std::min(firsthitfltlen,lasthitfltlen);
          TrkHelixUtils::findZFltlen(particle.traj(),_zent,entlen,0.1);
          _eventinfo.reco_p0=particle.momentum(entlen).mag();

          TrkHitVector const& hots = particle.hitVector();
          for(auto iter=hots.begin(); iter!=hots.end(); iter++)
          {
            const TrkHit *hitOnTrack = *iter;
            const mu2e::TrkStrawHit* trkStrawHit = dynamic_cast<const mu2e::TrkStrawHit*>(hitOnTrack);
            if(trkStrawHit)
            {
              //get the sim particles which caused this hit
              int stepPointVectorsIndex = trkStrawHit->index();
              const mu2e::PtrStepPointMCVector &stepPointVector = stepPointVectors->at(stepPointVectorsIndex);
              for(auto iterStepPoint=stepPointVector.begin(); iterStepPoint!=stepPointVector.end(); iterStepPoint++)
              {
                const art::Ptr<StepPointMC> stepPoint = *iterStepPoint;
                simParticles[stepPoint->simParticle()]++;
              }
            }
          }

          //find simparticle which matches the highest number of straw hits associated with the reconstructed track
          int largestNumberOfHits=0;
          art::Ptr<mu2e::SimParticle> matchingSimparticle;
          std::map<art::Ptr<mu2e::SimParticle>, int>::const_iterator simParticleIter;
          for(simParticleIter=simParticles.begin(); simParticleIter!=simParticles.end(); simParticleIter++)
          {
            if(simParticleIter->second>largestNumberOfHits)
            {
              largestNumberOfHits=simParticleIter->second;
              matchingSimparticle=simParticleIter->first;
            }
          }

          //access the simparticle found above
          if(largestNumberOfHits>0)
          {
            _eventinfo.simreco_endp=matchingSimparticle->endMomentum().vect().mag();

            //this sim particle may have been create in the previous simulation stage
            //therefore, a particle which has the creation code "primary", may not be a primary
            //but a continuation of a particle from a previous stage.
            //this particle from the previous stage is this particle's parent
            if(matchingSimparticle->creationCode()==ProcessCode::mu2ePrimary &&
               matchingSimparticle->hasParent()) matchingSimparticle=matchingSimparticle->parent();

            ProcessCode productionProcess = matchingSimparticle->creationCode();
            strncpy(_eventinfo.simreco_production_process, productionProcess.name().c_str(),99);

            int pdgId=matchingSimparticle->pdgId();
            strncpy(_eventinfo.simreco_particle, HepPID::particleName(pdgId).c_str(),99);

            _eventinfo.simreco_startp=matchingSimparticle->startMomentum().vect().mag();
            _eventinfo.simreco_startp_z=matchingSimparticle->startMomentum().vect().z();
            _eventinfo.simreco_pos[0]=matchingSimparticle->startPosition().x()-_detSysOrigin.x();
            _eventinfo.simreco_pos[1]=matchingSimparticle->startPosition().y()-_detSysOrigin.y();
            _eventinfo.simreco_pos[2]=matchingSimparticle->startPosition().z()-_detSysOrigin.z();

            std::string productionVolumeName="unknown volume";
            if(hasPhysicalVolumes && physicalVolumesMulti.product()!=NULL)
            {
              mu2e::PhysicalVolumeMultiHelper volumeMultiHelper(physicalVolumesMulti.product());
              productionVolumeName=volumeMultiHelper.startVolume(*matchingSimparticle).name();
            }
            strncpy(_eventinfo.simreco_production_volume, productionVolumeName.c_str(),99);
          }
        } //kalReps.size()>0
      } //kalReps.product()!=0
    } //event.getByLabel

    if(strlen(_eventinfo.simreco_particle)==0) //this mostly likely happened because the reconstructed track wasn't found
    {                                          //don't record such events
      std::cout<<"Will not record event "<<_eventinfo.event_number<<" of subrun "<<_eventinfo.subrun_number<<" in file "<<_eventinfo.filename<<", ";
      std::cout<<"since either the reco or the MC track is missing."<<std::endl;
      return;
    }

//check CRV veto
    GeomHandle<CosmicRayShield> CRS;

    art::Handle<CrvCoincidenceCheckResult> crvCoincidenceCheckResult;
    std::string crvCoincidenceInstanceName="";

    if(event.getByLabel(_crvCoincidenceModuleLabel,crvCoincidenceInstanceName,crvCoincidenceCheckResult))
    {
      const std::vector<CrvCoincidenceCheckResult::CoincidenceCombination> &coincidenceCombinations = crvCoincidenceCheckResult->GetCoincidenceCombinations();
      for(unsigned int i=0; i<coincidenceCombinations.size(); i++)
      {
        CRSScintillatorBarIndex barIndex=coincidenceCombinations[i]._counters[0];
        const CRSScintillatorBar &CRSbar = CRS->getBar(barIndex);
        int sectorNumber = CRSbar.id().getShieldNumber();
        int sectorType = CRS->getCRSScintillatorShield(sectorNumber).getSectorType();
        sectorType--;
        if(sectorType>=8 || sectorType<0) continue;

        _eventinfo.CRVveto_allSectors=true;
        _eventinfo.CRVveto[sectorType]=true;

        for(int j=0; j<3; j++)
        {
          double t=coincidenceCombinations[i]._time[j];
          barIndex=coincidenceCombinations[i]._counters[j];
          CLHEP::Hep3Vector pos = CRS->getBar(barIndex).getPosition()-_detSysOrigin;

          if(std::isnan(_eventinfo.firstCoincidenceHitTime) || t<_eventinfo.firstCoincidenceHitTime)
          {
            _eventinfo.firstCoincidenceHitTime=t;
            _eventinfo.firstCoincidenceHitSectorType=sectorType;
            for(int k=0; k<3; k++) _eventinfo.firstCoincidenceHitPos[k]=pos[k];
          }

          if(std::isnan(_eventinfo.CRVvetoTime[sectorType]) || t<_eventinfo.CRVvetoTime[sectorType])
          {
            _eventinfo.CRVvetoTime[sectorType]=t;
            for(int k=0; k<3; k++) _eventinfo.CRVvetoPos[sectorType][k]=pos[k];
          }
        }
      }
    }

//check for CRV step points
    std::vector<art::Handle<StepPointMCCollection> > CRVStepsVector;
    art::Selector selector(art::ProductInstanceNameSelector("CRV") &&
                           art::ProcessNameSelector("*"));

    event.getMany(selector, CRVStepsVector);
    for(size_t i=0; i<CRVStepsVector.size(); i++)
    {
      const art::Handle<StepPointMCCollection> &CRVSteps = CRVStepsVector[i];

      for(StepPointMCCollection::const_iterator iter=CRVSteps->begin(); iter!=CRVSteps->end(); iter++)
      {
        StepPointMC const& step(*iter);
        int PDGcode = step.simParticle()->pdgId();
//        if(abs(PDGcode)!=11 && abs(PDGcode)!=13) continue;
        if(abs(PDGcode)!=13) continue;   //ignore cases where a CRV hit is created by a particle other than a muon

        const CRSScintillatorBar &CRSbar = CRS->getBar(step.barIndex());

        int sectorNumber = CRSbar.id().getShieldNumber();
        int sectorType = CRS->getCRSScintillatorShield(sectorNumber).getSectorType();
        sectorType--;
        if(sectorType>=8 || sectorType<0) continue;
        // 0: R
        // 1: L
        // 2: T
        // 3: D
        // 4: U
        // 5,6,7: C1,C2,C3

        double t = step.time();
        CLHEP::Hep3Vector pos = step.position()-_detSysOrigin;
        CLHEP::Hep3Vector dir = step.momentum().unit();

        _eventinfo.CRVhit_allSectors =true;
        _eventinfo.CRVhit[sectorType]=true;
        if(_eventinfo.CRVhitTime[sectorType]>t || std::isnan(_eventinfo.CRVhitTime[sectorType]))
        {
          _eventinfo.CRVhitTime[sectorType]=t;
          for(int j=0; j<3; j++)
          {
            _eventinfo.CRVhitPos[sectorType][j]=pos[j];
            _eventinfo.CRVhitDir[sectorType][j]=dir[j];
          }
        }
      }
    }
*/

/***************/
/* only a test */
/*
    art::Handle<CrvRecoPulsesCollection> crvRecoPulsesCollection;
    event.getByLabel("CrvRecoPulses","",crvRecoPulsesCollection);
    for(CrvRecoPulsesCollection::const_iterator iter=crvRecoPulsesCollection->begin();
        iter!=crvRecoPulsesCollection->end(); iter++)
    {
      const CRSScintillatorBarIndex &barIndex = iter->first;
      const CRSScintillatorBar &CRSbar = CRS->getBar(barIndex);

      const CrvRecoPulses &crvRecoPulses = iter->second;
      for(int SiPM=0; SiPM<4; SiPM++)
      {
        const std::vector<CrvRecoPulses::CrvSingleRecoPulse> &pulseVector = crvRecoPulses.GetRecoPulses(SiPM);
        for(unsigned int i = 0; i<pulseVector.size(); i++)
        {
          const CrvRecoPulses::CrvSingleRecoPulse &pulse = pulseVector[i];
          int PEs=pulse._PEs;
          double time=pulse._leadingEdge;

          std::cout<<"run "<<event.id().run()<<"    subrun "<<event.id().subRun()<<"    event "<<event.id().event()<<"    barID "<<CRSbar.id()<<"    SiPM "<<SiPM<<"     PEs "<<PEs<<"      time "<<time<<std::endl;
        }
      }
    }
*/
/***************/

//FIXME
/*
    _tree->Fill();
*/
  }
}

using mu2e::CosmicAnalysis;
DEFINE_ART_MODULE(CosmicAnalysis)

