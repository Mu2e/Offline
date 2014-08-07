//
// A module to create CRV hits from StepPointMCs
//
// Based on the ideas found in 
// "A Monte Carlo simulation code for the
// development of WLS fiber readout plastic
// scintillator detector for the GRAPES-3 EAS
// array" by S.K. Gupta (gupta@grapes.tifr.res.in)
//
// $Id: MakeCRVHits_module.cc,v 1.1 2014/08/07 01:33:40 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:40 $
// 
// Original Author: Ralf Ehrlich

#include "CLHEP/Units/GlobalSystemOfUnits.h"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StepInstanceName.hh"
#include "MCDataProducts/inc/CRVHitCollection.hh"
#include "art/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h"

#include <string>

#include <TRandom3.h>
#include <TMath.h>

namespace mu2e 
{

  static const double c=2.99792458e2;  //in mm/ns
  static const double massE=0.510998910;  //in MeV/c^2
  static const double massMu=105.658367;  //in MeV/c^2

  class MakeCRVHits : public art::EDProducer 
  {

    public:
    explicit MakeCRVHits(fhicl::ParameterSet const& pset);
    void produce(art::Event& e);
    void beginJob();
    void endJob();

    private:
    std::string _g4ModuleLabel;
    double      _bluePhotonsPerMeV;
    double      _decayTimeScintillator;
    double      _decayTimeFiber;
    double      _indexOfRefractionScintillator;
    double      _indexOfRefractionFiber;
    double      _indexOfRefractionClad1;
    double      _indexOfRefractionClad2;
    double      _ETIRScintillator; //effective total internal reflectivity
    double      _ETIRFiber;        //effective total internal reflectivity
    double      _attenuationLengthScintillator;
    double      _attenuationLengthFiber;
    double      _titaniumOxideReflectivity;
    double      _fiberDistanceFromCenter;
    double      _fiberRadius;
    double      _photoDetectionEfficiency;
    TRandom3    _random3;
    double      _speedOfLightScintillator;
    double      _speedOfLightFiber;
    double      _criticalAngleScintillatorAir;
    double      _criticalAngleScintillatorClad2;
    double      _criticalAngleClad1Clad2;
    double      _criticalAngleFiberClad1;

    int         _longestAxis;       //that's the direction of the fiber
    int         _intermediateAxis;  //fibers sit left and right around the center of the intermediate axis
    int         _shortestAxis;      //fibers sit in the center of the shortest axis
  
    struct ScintillatorPhoton
    {
      CLHEP::Hep3Vector p0;
      CLHEP::Hep3Vector velocity;
      double t0;
    };
    struct FiberPhoton
    {
      CLHEP::Hep3Vector velocity;
      double l0, t0, r0;
      int fiberNumber;
    };
    struct PhotoElectron
    {
      double time;
      int fiberNumber;
      bool positiveSide;
    };

    CLHEP::Hep3Vector calculateVelocity(const CLHEP::Hep3Vector &momentum, const int &particleId);
    bool findBorderTime(const CLHEP::Hep3Vector &position0, const CLHEP::Hep3Vector &velocity, const double &time0, 
                        const CRSScintillatorBar &CRSbar, double &time1, double &time2);
    bool findBorderTime(const CLHEP::Hep3Vector &position0, const CLHEP::Hep3Vector &velocity, const double &time0, 
                        const CRSScintillatorBar &CRSbar, double &time1, double &time2, 
                        int &edgeNumber, bool &positiveSide);
    void propagateScintillatorPhotons(const std::vector<ScintillatorPhoton> &scintillatorPhotons,
                                      std::vector<FiberPhoton> &fiberPhotons,
                                      const CRSScintillatorBar &CRSbar);
    bool findReflection(CLHEP::Hep3Vector &velocity, const CRSScintillatorBar &CRSbar,
                        const int &edgeNumber, const bool &positiveSide);
    int checkFiberHits(const CLHEP::Hep3Vector &p0, const double &t0, 
                       const CLHEP::Hep3Vector &velocity, const CRSScintillatorBar &CRSbar, 
                       double &l, double &t2, int &fiberNumber,
                       CLHEP::Hep3Vector &reflectedVelocity);
    void propagateFiberPhotons(const std::vector<FiberPhoton> &scintillatorPhotons,
                               std::vector<PhotoElectron> &photoElectrons,
                               const CRSScintillatorBar &CRSbar);
  };

  MakeCRVHits::MakeCRVHits(fhicl::ParameterSet const& pset) :
    _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel","g4filter")),
    _bluePhotonsPerMeV(pset.get<double>("bluePhotonsPerMeV",7e3)), //7e3 //1.25e4 blue photons per MeV deposited energy
    _decayTimeScintillator(pset.get<double>("decayTimeScintillator",1.0)),                   //1ns
    _decayTimeFiber(pset.get<double>("decayTimeFiber",6.0)),                                 //6ns
    _indexOfRefractionScintillator(pset.get<double>("indexOfRefractionScintillator",1.59)),  //1.59
    _indexOfRefractionFiber(pset.get<double>("indexOfRefractionFiber",1.59)),                //1.59
    _indexOfRefractionClad1(pset.get<double>("indexOfRefractionClad1",1.49)),                //1.49
    _indexOfRefractionClad2(pset.get<double>("indexOfRefractionClad2",1.42)),                //1.42
    _ETIRScintillator(pset.get<double>("ETIRScintillator",0.93)),                            //0.93
    _ETIRFiber(pset.get<double>("ETIRFiber",0.9999)),                                        //0.9999
    _attenuationLengthScintillator(pset.get<double>("attenuationLengthScintillator",1000.0)),//1000mm
    _attenuationLengthFiber(pset.get<double>("attenuationLengthFiber",3500)),                //3500mm
    _titaniumOxideReflectivity(pset.get<double>("titaniumOxideReflectivity",0.90)),          //0.90
    _fiberDistanceFromCenter(pset.get<double>("fiberDistanceFromCenter",5.0)),               //5.0mm
    _fiberRadius(pset.get<double>("fiberRadius",0.658)),                                     //0.7mm-6%(clad thickness)
    _photoDetectionEfficiency(pset.get<double>("photoDetectionEfficiency",0.3))              //0.3
  {
    _speedOfLightScintillator = c/_indexOfRefractionScintillator;
    _speedOfLightFiber = c/_indexOfRefractionFiber;
    _criticalAngleScintillatorAir = asin(1/_indexOfRefractionScintillator);
    _criticalAngleScintillatorClad2 = asin(_indexOfRefractionClad2/_indexOfRefractionScintillator);
    _criticalAngleFiberClad1 = asin(_indexOfRefractionClad1/_indexOfRefractionFiber);
    _criticalAngleClad1Clad2 = asin(_indexOfRefractionClad2/_indexOfRefractionClad1);
    produces<CRVHitCollection>();
  }

  void MakeCRVHits::beginJob()
  {
  }

  void MakeCRVHits::endJob()
  {
  }

  void MakeCRVHits::produce(art::Event& event) 
  {
    std::unique_ptr<CRVHitCollection> crvHits(new CRVHitCollection);

    GeomHandle<CosmicRayShield> CRS;
    StepInstanceName CRVInstance(StepInstanceName::CRV);
    art::Handle<StepPointMCCollection> CRVSteps;
    event.getByLabel(_g4ModuleLabel,CRVInstance.name(),CRVSteps);

    for(StepPointMCCollection::const_iterator iter=CRVSteps->begin(); iter!=CRVSteps->end(); iter++)
    {
      StepPointMC const& step(*iter);
      cet::map_vector_key trackId=step.simParticle()->id();
      int particleId=step.simParticle()->pdgId();
      if(abs(particleId)!=11 && abs(particleId)!=13) continue;

      const CRSScintillatorBar &CRSbar = CRS->getBar(step.barIndex());

      //find the muon/electron track segment (enter and exit times in the counter)
      CLHEP::Hep3Vector velocity = calculateVelocity(step.momentum(),particleId);
      double time0=step.time();
      double time1, time2;
      if(step.momentum().mag()!=0)
      {
        if(!findBorderTime(step.position(), velocity, time0, CRSbar, time1, time2)) 
        {
          std::cout<<"Looks like the StepPointMC point was outside of the scintillator bar. ";
          std::cout<<"This shouldn't happen."<<std::endl; 
          continue;
        }
      }
      else
      {
        time1=time0;
        time2=time0;
      }

      //modify time1, if there was a previous step in this scintillator bar
      if(iter!=CRVSteps->begin())
      {
        StepPointMCCollection::const_iterator iterPrevStep = iter;
        iterPrevStep--;
        if(iterPrevStep->barIndex()==step.barIndex()
           && trackId==iterPrevStep->simParticle()->id())
        {
          double timePrev = iterPrevStep->time();
          double newTime1 = (time0+timePrev)/2.0;
          if(newTime1>time1) time1=newTime1;
        }
      }

      //modify time2, if there is a following step in this scintillator bar
      StepPointMCCollection::const_iterator iterNextStep = iter;
      iterNextStep++;
      if(iterNextStep!=CRVSteps->end())
      {
        if(iterNextStep->barIndex()==step.barIndex()
           && trackId==iterNextStep->simParticle()->id())
        {
          double timeNext = iterNextStep->time();
          double newTime2 = (time0+timeNext)/2.0;
          if(newTime2<time2) time2=newTime2;
        }
      }

      const std::vector<double> &barHalfLengths = CRSbar.getHalfLengths();
      //fibers go along the longest axis
      _longestAxis=std::distance(barHalfLengths.begin(),std::max_element(barHalfLengths.begin(),barHalfLengths.end()));
      //fibers are located at the center of the shortest axis
      _shortestAxis=std::distance(barHalfLengths.begin(),std::min_element(barHalfLengths.begin(),barHalfLengths.end()));
      //the two fibers are located "left" and "right" around the center of the intermediate axis
      for(_intermediateAxis=0; _intermediateAxis<3; _intermediateAxis++)
      {
        if(_intermediateAxis!=_shortestAxis && _intermediateAxis!=_longestAxis) break;
      }

      //the track section which will be considered now is between time1 and time2 around the current step point
      //blue photons will be emitted along this tracks in equal distances
      int numberOfPhotons=static_cast<int>(_bluePhotonsPerMeV*step.totalEDep());
      if(numberOfPhotons==0) continue;

      std::vector<ScintillatorPhoton> scintillatorPhotons;
      std::vector<FiberPhoton> fiberPhotons;
      std::vector<PhotoElectron> photoElectrons;
      scintillatorPhotons.reserve(numberOfPhotons);
      fiberPhotons.reserve(numberOfPhotons);
      photoElectrons.reserve(numberOfPhotons);

      double timeInterval=(time2-time1)/numberOfPhotons;
      for(int i=0; i<numberOfPhotons; i++)
      {
        double t0 = time1+i*timeInterval+0.5*timeInterval;                  //time when the energy gets deposited  
        CLHEP::Hep3Vector p0 = step.position() + velocity*(t0-time0);       //position where the energy gets deposited
        double r = _random3.Rndm();
        if(r==1.0) r=0.0;
        t0 += -_decayTimeScintillator*log(1.0-r);             //time when the scintillator emits the blue photon
        //uniform emission direction, so the setup of the Mu2e coordinate system is not relevant
        double theta = acos(1.0-2.0*_random3.Rndm());         //theta of emitted photon
        double phi   = 2.0*TMath::Pi()*_random3.Rndm();       //phi of the emitted photon
        CLHEP::Hep3Vector velocityPhoton;
        velocityPhoton.setRThetaPhi(_speedOfLightScintillator, theta, phi);

        ScintillatorPhoton p;
        p.t0 = t0;
        p.p0 = p0;
        p.velocity = velocityPhoton;
        scintillatorPhotons.push_back(p);
      }

      propagateScintillatorPhotons(scintillatorPhotons,fiberPhotons,CRSbar);
      propagateFiberPhotons(fiberPhotons,photoElectrons,CRSbar);

      std::set<CRVHit> &crvHitSet = (*crvHits)[step.barIndex()];
      std::vector<PhotoElectron>::const_iterator PEiter;
      for(PEiter=photoElectrons.begin(); PEiter!=photoElectrons.end(); PEiter++)
      {
        CRVHit hit(PEiter->time, PEiter->fiberNumber, PEiter->positiveSide?1:0);
        crvHitSet.insert(hit);
      }
    }

    event.put(std::move(crvHits));
  } // end produce

  CLHEP::Hep3Vector MakeCRVHits::calculateVelocity(const CLHEP::Hep3Vector &momentum, const int &particleId)
  {
    double mass = abs(particleId)==13?massMu:massE;           //mass in MeV/c^2
    double velocity = 1/(pow(mass/momentum.mag(),2)+1);       //velocity in units of c
    velocity *= c;                                            //velocity in mm/ns
    CLHEP::Hep3Vector velocityVector = momentum.unit()*velocity;
    return velocityVector;
  }

  bool MakeCRVHits::findBorderTime(const CLHEP::Hep3Vector &position0, const CLHEP::Hep3Vector &velocity, 
                                   const double &time0, 
                                   const CRSScintillatorBar &CRSbar, double &time1, double &time2)
  {
    int dummyEdgeNumber;  //not needed if we simply want to find the time when the track hits the border
    bool dummyPositiveSide;
    return findBorderTime(position0, velocity, time0, CRSbar, time1, time2, dummyEdgeNumber, dummyPositiveSide);
  }

  bool MakeCRVHits::findBorderTime(const CLHEP::Hep3Vector &position0, const CLHEP::Hep3Vector &velocity, 
                                   const double &time0, 
                                   const CRSScintillatorBar &CRSbar, double &time1, double &time2,
                                   int &edgeNumber, bool &positiveSide)
  {
    const CLHEP::Hep3Vector &barPosition = CRSbar.getPosition();
    const std::vector<double> &barHalfLengths = CRSbar.getHalfLengths();
    double deltaTime1[3], deltaTime2[3];
    bool positiveSideArray[3]={true,true,true};
    for(int i=0; i<3; i++)  //all 3 dimensions
    {
      double p0 = position0[i];
      double v  = velocity[i];
      double p1 = barPosition[i]-barHalfLengths[i];
      double p2 = barPosition[i]+barHalfLengths[i];
      if(v!=0)
      {
        //p1 = p0 + v*deltaTime1 
        //p2 = p0 + v*deltaTime2
        deltaTime1[i]=(p1-p0)/v;
        deltaTime2[i]=(p2-p0)/v;
        if(deltaTime1[i]>deltaTime2[i])
        {
          std::swap(deltaTime1[i],deltaTime2[i]); //keep all times in ascending order
          positiveSideArray[i]=false;
        }
      }
      else  //no velocity component in the i-th dimension, 
            //which means that the track will never cross the borders in the i-th dimension
      {
        if(p0<p1 || p0>p2) return(false);  //track outside of scintillator bar. this shouldn't happen.
        deltaTime1[i]=-INFINITY;
        deltaTime2[i]=+INFINITY;
      }
    }

    //find overlaps of the passage times for all three dimensions
    //if there is an overlap, the trajectory went through the scintillator bar
    //and the entrance and exit time can be determined
    //example:
    //T1:     +++++++++++
    //T2:          +++++++++++++
    //T3:       +++++++
    //overlap:     ====
    //----------------------------------------------->t
    
    time1 = *std::max_element(deltaTime1,deltaTime1+3) + time0;
    time2 = *std::min_element(deltaTime2,deltaTime2+3) + time0;
    if(time1>time2) return(false); //no overlap found (track was not in scintillator bar). this shouldn't happen.
    if(time1>time0 && time1-time0<1e-9) time1=time0; //small corrections due to rounding errors which may occur
    if(time2<time0 && time0-time2<1e-9) time2=time0; //small corrections due to rounding errors which may occur
    if(time1>time0 || time2<time0) return(false);    //start point was not in scintillator bar. this shouldn't happen.
    if(isinf(time1) || isinf(time2)) return(false);  //this shouldn't happen.

    edgeNumber = std::min_element(deltaTime2,deltaTime2+3)-deltaTime2;   //which entry in deltaTime2 was the smallest
                                                                         //i.e. the one responsible for a potential 
                                                                         //reflection (= the index of the bar wall).
    positiveSide = positiveSideArray[edgeNumber];
    return(true);
  }

  bool MakeCRVHits::findReflection(CLHEP::Hep3Vector &velocity, const CRSScintillatorBar &CRSbar,
                                   const int &edgeNumber, const bool &positiveSide)
  {
    if(_longestAxis==edgeNumber) return(false); //no internal reflection on the plastic end caps

    CLHEP::Hep3Vector velocityBefore=velocity;
    velocity[edgeNumber]=-velocity[edgeNumber];
    /*
    the negative of the "velocityBefore" vector has to be taken when calculating the angle between the
    incoming and reflected light, so that both vectors are pointing away from the surface
    */
    double incidentAngle = velocity.angle(-velocityBefore)/2.0;
    bool photonEscapes=false;
    if(incidentAngle<=_criticalAngleScintillatorAir) photonEscapes=true;  //no internal reflection; photon may escape
    else
    {
      if(_random3.Rndm()>_ETIRScintillator) 
      photonEscapes=true;   //no internal reflection due to rough surface; photon may escape
    }
    if(photonEscapes)       //escaped photons may be reflected on the titanium oxide layer
    {
      if(_random3.Rndm()>_titaniumOxideReflectivity) return(false); //no reflection on titanium oxide 
                                                                    //due to its reflectivity 
      //direction from the point of view of the titanium oxide surface (with z component point away from the surface)
      double theta=asin(_random3.Rndm());
      double phi=2.0*TMath::Pi()*_random3.Rndm();
      velocity.setRThetaPhi(_speedOfLightScintillator, theta, phi);
      //need to swap components depending on which border the light was reflected.
      double v[3]={velocity.x(),velocity.y(),velocity.z()};
      //the z component (=normal component) from the titanium oxide point of view (0.0 ... +_speedOfLightScintillator) 
      //is always positive. needs to be swapped with the edgeNumber-th component.
      //the other two components can be defined arbitrarily, since phi is distributed uniformly.
      std::swap(v[edgeNumber],v[2]);
      //the direction of this component needs to be changed if the reflection happens on the positive end of the bar
      if(positiveSide) v[edgeNumber]*=-1.0;
      velocity.set(v[0],v[1],v[2]);
    }
    return(true);
  }
  
  void MakeCRVHits::propagateScintillatorPhotons(const std::vector<ScintillatorPhoton> &scintillatorPhotons,
                                                 std::vector<FiberPhoton> &fiberPhotons,
                                                 const CRSScintillatorBar &CRSbar)
  {
    std::vector<ScintillatorPhoton>::const_iterator iter;
    for(iter=scintillatorPhotons.begin(); iter!=scintillatorPhotons.end(); iter++)
    {
      CLHEP::Hep3Vector velocity = iter->velocity;
      CLHEP::Hep3Vector p0 = iter->p0;
      double t0 = iter->t0;

      double r=_random3.Rndm();
      if(r==1.0) r=0.0;
      double survivalLength=-_attenuationLengthScintillator*log(1-r);
      double photonPathLength=0;
      for(int iReflection=0; iReflection<150; iReflection++)  //up to 150 reflections
      {
        double lFiber,t2;
        int fiberNumber; 
        CLHEP::Hep3Vector reflectedVelocity;
        int fiberHit=checkFiberHits(p0,t0,velocity,CRSbar,lFiber,t2,fiberNumber, reflectedVelocity);
        if(fiberHit==2) //photon got absorbed by fiber
        {
          CLHEP::Hep3Vector p2 = p0 + velocity*(t2-t0);   //p2, t2 are location and time where the photon hits
                                                          //the fiber (origin of the green photon).
          photonPathLength+=(p2-p0).mag();
          if(photonPathLength>survivalLength) break;      //photon lost due to attenuation

          t2 += -_decayTimeFiber*log(1.0-r);              //time when the fiber emits the green photon
          //uniform emission direction of the green photon
          double theta = acos(1.0-2.0*_random3.Rndm());   //theta of emitted green photon
          double phi   = 2.0*TMath::Pi()*_random3.Rndm(); //phi of the emitted green photon
          CLHEP::Hep3Vector velocityFiber;
          velocityFiber.setRThetaPhi(_speedOfLightFiber, theta, phi);

          //the position along the photon path inside of the fiber 
          //where the blue photon gets absorbed and the green photon gets emitted is not known.
          //l0 is assumed to be where the photon first hits the surface of the fiber 
          //(the error will be small due to the small fiber radius).
          //r0 is assumed to uniformily distributed between 0 and fiberRadius.
          //this is not true, but i don't have a better method for the time being
          double r0    = _fiberRadius*_random3.Rndm();

          //put everything into a new vector of photons absorbed by fiber
          FiberPhoton p;
          p.t0 = t2;
          p.l0 = lFiber;
          p.r0 = r0;
          p.velocity = velocityFiber;
          p.fiberNumber = fiberNumber;
          fiberPhotons.push_back(p);
          break;   //end propagation, since the photon has been absorbed
        }
        if(fiberHit==1) //photon got reflected on fiber
        {
          CLHEP::Hep3Vector p2 = p0 + velocity*(t2-t0); //p2, t2 are location and time where the photon hits the fiber
          photonPathLength+=(p2-p0).mag();
          if(photonPathLength>survivalLength) break;    //photon lost due to attenuation
          velocity=reflectedVelocity;
          t0=t2;
          p0=p2;
          continue;  //restart the loop to find the next reflection
        }

        //photon has not been absorbed by fiber, need to calculate next reflection
        double t1;    //t1 is dummy variable
        int edgeNumber;
        bool positiveSide;
        if(!findBorderTime(p0, velocity, t0, CRSbar, t1, t2, edgeNumber, positiveSide)) break; //shouldn't happen
        CLHEP::Hep3Vector p2 = p0 + velocity*(t2-t0);   //p2, t2 are location and time where the photon hits
                                                        //the edge of the scintillator bar.

        photonPathLength+=(p2-p0).mag();
        if(photonPathLength>survivalLength) break;      //photon lost due to attenuation

        if(!findReflection(velocity, CRSbar, edgeNumber, positiveSide)) break; //photon escapes
        t0=t2;
        p0=p2;
        //restart the loop to find the next reflection
      }
    }
  }

  int MakeCRVHits::checkFiberHits(const CLHEP::Hep3Vector &p0, const double &t0, 
                                  const CLHEP::Hep3Vector &velocity, const CRSScintillatorBar &CRSbar, 
                                  double &l, double &t2, int &fiberNumber, 
                                  CLHEP::Hep3Vector &reflectedVelocity)
  {
    const CLHEP::Hep3Vector &barPosition = CRSbar.getPosition();
    const std::vector<double> &barHalfLengths = CRSbar.getHalfLengths();
    
    CLHEP::Hep3Vector vectorLongestAxis, vectorIntermediateAxis;
    vectorLongestAxis[_longestAxis]=barHalfLengths[_longestAxis];
    vectorIntermediateAxis[_intermediateAxis]=_fiberDistanceFromCenter;

    //start point of the two fibers
    CLHEP::Hep3Vector fiberP0[2];
    fiberP0[0]=barPosition-vectorLongestAxis-vectorIntermediateAxis;
    fiberP0[1]=barPosition-vectorLongestAxis+vectorIntermediateAxis;

    bool capturedByFiber=false;
    bool reflectedByFiber=false;
    double deltaT=NAN;
    for(int i=0; i<2; i++)   //loop over both fibers
    {
      /*
      the equation of the photon track (split into components,
      where the actual direction of the 1st, 2nd, 3rd component is determined by
      the orientation of the scintillator bar)
      (1) trackP_1 = trackP0_1 + trackVelocity_1 * deltaT [_shortestAxis]
      (2) trackP_2 = trackP0_2 + trackVelocity_2 * deltaT [_intermediateAxis]
      (3) trackP_3 = trackP0_3 + trackVelocity_3 * deltaT [_longestAxis]

      the equation of the cylinder representing the fiber (split into components, 
      where the actual direction of the 1st, 2nd, 3rd component is determined by
      the orientation of the scintillator bar)
      (4) fiberP_1 = fiberP0_1 + R*sin(phi)  [_shortestAxis]
      (5) fiberP_2 = fiberP0_2 + R*cos(phi)  [_intermediateAxis]
      (6) fiberP_3 = fiberP0_3 + l           [_longestAxis]
      */

      //start with (1)=(4) and solve for deltaT  --> equation (7)
      //use (2)=(5), plug in deltaT from above and solve for phi
      double V=velocity[_intermediateAxis]/velocity[_shortestAxis];
      double D=((fiberP0[i][_intermediateAxis]-p0[_intermediateAxis])-V*(fiberP0[i][_shortestAxis]-p0[_shortestAxis]))
               /(_fiberRadius*sqrt(V*V+1));
      if(fabs(D)>1.0) continue;    //no solution for phi. photon track doesn't interect fiber cylinder.
      //phi has generally two solutions
      double phi1, phi2;
      phi1=asin(D)-atan2(-1.0,V);
      if(D>=0.0) phi2=(TMath::Pi()-asin(D))-atan2(-1.0,V);
      else phi2=(-TMath::Pi()-asin(D))-atan2(-1.0,V);

      //plug both solutions of phi into (7) to find two solutions for deltaT
      double deltaT1=((fiberP0[i][_shortestAxis]-p0[_shortestAxis])+_fiberRadius*sin(phi1))/velocity[_shortestAxis];
      double deltaT2=((fiberP0[i][_shortestAxis]-p0[_shortestAxis])+_fiberRadius*sin(phi2))/velocity[_shortestAxis];
      if(deltaT1<=0 || deltaT2<=0) continue;  //this shouldn't happen, since it would mean
                                              //that the photon came from the direction of the fiber

      //use the smallest deltaT (and associated phi), since this gives the time when the photon first hits the fiber
      double deltaT_tmp=deltaT1<deltaT2?deltaT1:deltaT2;
      double phi_tmp=deltaT1<deltaT2?phi1:phi2;

      //use deltaT in (3) to find the 3rd component (=longest axis) of the hit position
      double pLongestAxis=p0[_longestAxis]+velocity[_longestAxis]*deltaT_tmp;
      //use (3)=(6) with the result from above to find the position along the fiber 
      //relative to the starting position of the fiber
      double l_tmp=pLongestAxis-fiberP0[i][_longestAxis];
      //check whether its within the fiber length
      if(l_tmp<=0 || l_tmp>=2.0*barHalfLengths[_longestAxis]) continue; //photon track intersected fiber outside of bar
      //the photon has hit the fiber

      //check whether it gets reflected on the fiber clad
      /*
      the normal vector on the surface should be given by the difference between 
      the impact point on the surface (which is fiberP from equations (4)-(6))
      and the point on the center of the fiber "underneath" the impact point
      (which would be fiberP from equations (4)-(6) with R=0) 
      (4) fiberP_1 = fiberP0_1 + R*sin(phi)  [_shortestAxis]
      (5) fiberP_2 = fiberP0_2 + R*cos(phi)  [_intermediateAxis]
      (6) fiberP_3 = fiberP0_3 + l           [_longestAxis]
      */
      CLHEP::Hep3Vector fiberSurfaceNormalVector;
      fiberSurfaceNormalVector[_shortestAxis]=_fiberRadius*sin(phi_tmp);
      fiberSurfaceNormalVector[_intermediateAxis]=_fiberRadius*cos(phi_tmp);
      fiberSurfaceNormalVector[_longestAxis]=0;

      /*
      the negative of the velocity vector has to be taken when calculating the angle between the
      incoming light and normal vector, so that both vectors are pointing away from the surface
      */
      double incidentAngle = fiberSurfaceNormalVector.angle(-velocity);
      if(incidentAngle>_criticalAngleScintillatorClad2)
      {
        //photon gets reflected on the fiber clad
        //projection of the negative of the velocity vector on the normal vector
        CLHEP::Hep3Vector velocityProjected = fiberSurfaceNormalVector.dot(-velocity)*fiberSurfaceNormalVector;
        //reflectedVelocity = 2*velocityProjected - (-velocity)
        reflectedVelocity = 2.0*velocityProjected + velocity;
        if(deltaT>deltaT_tmp || isnan(deltaT))
        {
          reflectedByFiber=true;
          capturedByFiber=false;
          deltaT=deltaT_tmp;
          t2=t0+deltaT;
          l=l_tmp;
          fiberNumber=i;
        }
      }
      else
      {      
        //photon gets captured
        if(deltaT>deltaT_tmp || isnan(deltaT))
        {
          reflectedByFiber=false;
          capturedByFiber=true;
          deltaT=deltaT_tmp;
          t2=t0+deltaT;
          l=l_tmp;
          fiberNumber=i;
        }
      }
    }

    if(reflectedByFiber) return(1);
    if(capturedByFiber) return(2);
    return 0;
  } 

  void MakeCRVHits::propagateFiberPhotons(const std::vector<FiberPhoton> &fiberPhotons,
                                          std::vector<PhotoElectron> &photoElectrons,
                                          const CRSScintillatorBar &CRSbar)
  {
    const std::vector<double> &barHalfLengths = CRSbar.getHalfLengths();

    std::vector<FiberPhoton>::const_iterator iter;
    for(iter=fiberPhotons.begin(); iter!=fiberPhotons.end(); iter++)
    {
      double r=_random3.Rndm();
      if(r==1.0) r=0.0;
      double survivalLength=-_attenuationLengthFiber*log(1-r);
      double photonPathLength=0;

      const CLHEP::Hep3Vector &velocity = iter->velocity;
      double t0 = iter->t0;
      double l0 = iter->l0;
      double r0 = iter->r0;
      int fiberNumber = iter->fiberNumber;

      //first step is the track segment from r=r0 to r=fiberRadius

      //theta of the velocity vector is the azimuth angle with 0 being in positive direction along the fiber
      //when the photon hits the fiber border at r=fiberRadius after coming from r=r0,
      //the azimuth angle of the incoming track (with respect to the normal of the fiber border) 
      //as seen from the impact point is thetaFiber
      double thetaImpact=fabs(TMath::PiOver2()-velocity.theta());

      //phi of the velocity vector is the polar angle in the plane of the fiber cross section.
      //when the photon hits the fiber border at r=fiberRadius after coming from r=r0,
      //the polar angle of the incoming track (with respect to the normal of the fiber border) 
      //as seen from the impact point is phiFiber (calculated from the law of sines)
      double phiImpact=asin(sin(velocity.phi())*r0/_fiberRadius);

      //distance between the point where the photon was emitted and the first reflection point
      double thirdAngle=TMath::Pi()-(velocity.phi()+phiImpact);
      double s1=sqrt(r0*r0+_fiberRadius*_fiberRadius-2.0*r0*_fiberRadius*cos(thirdAngle)); //distance in the
                                                                                           //plane of the fiber
                                                                                           //cross section
      double s2=s1/tan(velocity.theta());    //distance along the fiber
      bool positiveDirection;

      //special cases where the photon's track hits the SiPM right away without reflection
      if(l0+s2<=0)  //includes +/-inf
      {
        //going in negative direction
        positiveDirection=false;
        double firstDistance=sqrt(s1*s1+l0*l0);
        photonPathLength+=firstDistance;
        if(photonPathLength>survivalLength) continue;   //photon lost due to attenuation
      }
      else if(l0+s2>=2.0*barHalfLengths[_longestAxis])  //includes +/-inf
      {
        //going in positive direction
        positiveDirection=true;
        double firstDistance=sqrt(s1*s1+pow(2.0*barHalfLengths[_longestAxis]-l0,2));
        photonPathLength+=firstDistance;
        if(photonPathLength>survivalLength) continue;   //photon lost due to attenuation
      }
      else
      {
        //normal case: the first track segment goes to the border of the fiber to get reflected
        double firstDistance=sqrt(s1*s1+s2*s2);
        photonPathLength+=firstDistance;
        if(photonPathLength>survivalLength) continue;  //photon lost due to attenuation

        CLHEP::Hep3Vector velocityImpact, normalVectorImpact;
        velocityImpact.setRThetaPhi(_speedOfLightFiber, thetaImpact, phiImpact);
        normalVectorImpact.setRThetaPhi(_speedOfLightFiber, 0.0, 0.0);

        double impactAngle = normalVectorImpact.angle(velocityImpact);
        if(impactAngle<_criticalAngleFiberClad1)
        { 
          //photon escapes fiber, since the angle is too small for internal reflection.
          //photon is now in inner clad (clad1)
          double impactAngleClad=asin(_indexOfRefractionFiber/_indexOfRefractionClad1*sin(impactAngle));
          if(impactAngleClad<_criticalAngleClad1Clad2) continue; //photon escapes clad1, since the angle is too 
                                                                 //small for internal reflection.
          //the additional travel distance of the photons in the fiber clad will be ignored
          //the error will be small due to the small thickness of the clad
        }
        //distance along the fiber covered between two reflections
        double thirdAngleSegment=TMath::Pi()-2.0*phiImpact;
        double s1Segment=_fiberRadius
                        *sqrt(2.0-2.0*r0*_fiberRadius*cos(thirdAngleSegment)); //distance in the plane of the fiber
        double s2Segment=s1Segment/tan(velocity.theta());                      //distance along the fiber
        double segmentDistance=sqrt(s1Segment*s1Segment+s2Segment*s2Segment); 

        //number of single segments needed to reach a SiPM
        if(s2Segment==0) continue;  //this photon will never reach a SiPM
        double requiredDistanceAlongFiber;
        if(s2Segment<0) requiredDistanceAlongFiber=l0-fabs(s2);
        else requiredDistanceAlongFiber=(2.0*barHalfLengths[_longestAxis]-l0)-fabs(s2);

        positiveDirection=(s2Segment>0?true:false);

        //total tracklength after it reached the SiPM
        double numberOfSegments=requiredDistanceAlongFiber/fabs(s2Segment);
        photonPathLength+=numberOfSegments*segmentDistance;;
        if(photonPathLength>survivalLength) continue;  //photon lost due to attenuation
      
        //number of reflections needed to reach a SiPM
        int numberOfReflections=ceil(numberOfSegments);
        //probability of internal reflection is less than 1.
        //total probability for n reflections is p^n.
        double totalReflectionProbability=pow(_ETIRFiber,numberOfReflections);
        if(_random3.Rndm()>totalReflectionProbability) continue; //photon escaped due to to imperfect internal reflection
      }

      //the photon reached the SiPM
      double arrivalTime=t0 + photonPathLength/_speedOfLightFiber;

      if(_random3.Rndm()>_photoDetectionEfficiency) continue; //photon could not produce a photo electron
      
      PhotoElectron p;
      p.time=arrivalTime;
      p.fiberNumber=fiberNumber;
      p.positiveSide=positiveDirection;
      photoElectrons.push_back(p);
    }
  }

} // end namespace mu2e

using mu2e::MakeCRVHits;
DEFINE_ART_MODULE(MakeCRVHits)
