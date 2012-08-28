// Pixel digitization: create ExtMonFNALRawHits and associated truth.
//
// $Id: ExtMonFNALHitMaker_module.cc,v 1.2 2012/08/28 05:05:52 gandr Exp $
// $Author: gandr $
// $Date: 2012/08/28 05:05:52 $
//
// Original author Andrei Gaponenko
//

#include <string>
#include <cmath>
#include <memory>
#include <queue>
#include <iostream>

#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "cetlib/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/Ptr.h"

#include "RecoDataProducts/inc/ExtMonFNALRawHit.hh"
#include "RecoDataProducts/inc/ExtMonFNALRawHitCollection.hh"
#include "MCDataProducts/inc/ExtMonFNALHitTruthAssn.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/ExtMonFNALSimHit.hh"
#include "MCDataProducts/inc/ExtMonFNALSimHitCollection.hh"

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALSensor.hh"
#include "ExtinctionMonitorFNAL/Digitization/inc/SiliconProperties.hh"
#include "ExtinctionMonitorFNAL/Digitization/inc/PixelToTCircuit.hh"
#include "ExtinctionMonitorFNAL/Digitization/inc/PixelNoise.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/ExtMonFNALConditions.hh"
#include "SeedService/inc/SeedService.hh"

namespace mu2e {

  //================================================================
  struct PixelTimedChargeDeposit {
    double time;
    double charge;
    art::Ptr<SimParticle> particle;

    PixelTimedChargeDeposit(double t, double c, const art::Ptr<SimParticle>& p)
      : time(t), charge(c), particle(p)
    {}

    // We accumulated deposits in a priority_queu, and want earlier times
    // to come out first.  Thus the inverted less-than definition:
    bool operator<(const PixelTimedChargeDeposit& b) const {
      return b.time < this->time;
    }
  };

  std::ostream& operator<<(std::ostream& os, const PixelTimedChargeDeposit& dep) {
    return os<<"PTD("<<dep.time<<", "<<dep.charge<<", "<<dep.particle<<")";
  }

  // Queue ordered by time
  class PixelChargeHistory : public std::priority_queue<PixelTimedChargeDeposit> {};

  class PixelChargeCollection : public std::map<ExtMonFNALPixelId,PixelChargeHistory> {};

  //================================================================
  class ExtMonFNALHitMaker : public art::EDProducer {

  public:
    explicit ExtMonFNALHitMaker(fhicl::ParameterSet const& pset)
      : inputModuleLabel_(pset.get<std::string>("inputModuleLabel"))
      , inputInstanceName_(pset.get<std::string>("inputInstanceName"))
      , nclusters_(pset.get<unsigned>("numClustersPerHit"))
      , maxToT_(pset.get<unsigned>("maxToT"))
      , discriminatorThreshold_(pset.get<double>("discriminatorThreshold"))
      , qCalib_(pset.get<double>("qCalib"))
      , totCalib_(pset.get<int>("totCalib"))

      , eng_(createEngine(art::ServiceHandle<SeedService>()->getSeed()))

      , gaussian_(eng_)

      , noise_(eng_,
               pset.get<double>("pixelNoisePerBC"),
               pset.get<int>("noiseClockMin"),
               pset.get<int>("noiseClockMax"))

      , extMon_(0)
      , cond_(0)
    {
      produces<ExtMonFNALRawHitCollection>();
      produces<ExtMonFNALHitTruthAssn>();

      if(nclusters_ < 2) {
        throw cet::exception("CONFIG")
          <<"ExtMonFNALHitMaker: the numClustersPerHit FHICL parameter should be greater than 1\n";
      }
    }

    virtual void produce(art::Event& evt);
    virtual void beginRun(art::Run& run);

  private:
    std::string inputModuleLabel_;
    std::string inputInstanceName_;

    // The number of charge clusters per one simhit
    unsigned nclusters_;

    // Dynamic range of readout ToT
    int maxToT_;

    double discriminatorThreshold_;
    double qCalib_;
    int    totCalib_;

    art::RandomNumberGenerator::base_engine_t& eng_;
    CLHEP::RandGaussQ gaussian_;

    SiliconProperties siProps_;

    PixelNoise noise_;

    // Non-owning pointers to the geometry and conditions objects. The
    // current Mu2e infrastructure does not allow the use of a Handle
    // as a class member.
    const ExtMonFNAL::ExtMon *extMon_;
    const ExtMonFNALConditions *cond_;

    void collectIonization(PixelChargeCollection *pixcharges,
                           const ExtMonFNALSimHitCollection& simhits);

    void collectIonization(PixelChargeCollection *pixcharges,
                           const ExtMonFNALSimHit& hit);

    void addCharge(PixelChargeCollection *pixcharges,
                   ExtMonFNALSensorId sid,
                   double time,
                   double charge,
                   double x_ro,
                   double y_ro,
                   const art::Ptr<SimParticle>& particle);

    // The input pixcharges collection gets eaten by this call
    void discriminate(ExtMonFNALRawHitCollection *outhits,
                      ExtMonFNALHitTruthAssn *outtruth,
                      art::ProductID outhitsPID,
                      const art::EDProductGetter* outhitsGetter,
                      PixelChargeCollection& pixcharges);

    // The input PixelChargeHistory gets eaten by this call
    void discriminate(ExtMonFNALRawHitCollection *outhits,
                      ExtMonFNALHitTruthAssn *outtruth,
                      art::ProductID outhitsPID,
                      const art::EDProductGetter* outhitsGetter,
                      const ExtMonFNALPixelId& pix,
                      PixelChargeHistory& ch);

    int timeStamp(double time) const {
      return (time - cond_->t0())/cond_->clockTick();
    }

  };

  //================================================================
  void ExtMonFNALHitMaker::beginRun(art::Run& run) {
    GeomHandle<ExtMonFNAL::ExtMon> emf;
    extMon_ = &*emf;

    ConditionsHandle<ExtMonFNALConditions> cond("ignored");
    cond_ = &*cond;

    const double sensorThickness = 2*emf->sensor().halfSize()[2];

    const double electricField = std::abs(cond->biasVoltage())/sensorThickness;

    siProps_.setConditions(cond->temperature(), electricField);

    std::cout<<"ExtMonFNALHitMaker:"
             <<"  electronHolePairsPerEnergy = "<<siProps_.electronHolePairsPerEnergy()
             <<", electronDriftMobility = "<<siProps_.electronDriftMobility()
             <<", electronDiffusionConstant = "<<siProps_.electronDiffusionConstant()
             <<std::endl;
  }

  //================================================================
  void ExtMonFNALHitMaker::produce(art::Event& event) {
    std::auto_ptr<ExtMonFNALRawHitCollection> outHits(new ExtMonFNALRawHitCollection());
    std::auto_ptr<ExtMonFNALHitTruthAssn> outTruth(new ExtMonFNALHitTruthAssn());

    art::Handle<ExtMonFNALSimHitCollection> ih;
    event.getByLabel(inputModuleLabel_, inputInstanceName_, ih);
    const ExtMonFNALSimHitCollection& simhits(*ih);

    // N.B.: to reduce memory footprint we can digitize one readout chip at a time
    // instead of doing the whole detector in one go.  Would be an extra loop here.
    PixelChargeCollection pixcharges;
    collectIonization(&pixcharges, simhits);

    const art::ProductID hitsPID = getProductID<ExtMonFNALRawHitCollection>(event);
    const art::EDProductGetter *hitsGetter = event.productGetter(hitsPID);
    discriminate(&*outHits, &*outTruth, hitsPID, hitsGetter, pixcharges);

    noise_.add(&*outHits);

    event.put(outHits);
    event.put(outTruth);
  }

  //================================================================
  void ExtMonFNALHitMaker::collectIonization(PixelChargeCollection *pixcharges,
                                             const ExtMonFNALSimHitCollection& simhits)
  {
    for(ExtMonFNALSimHitCollection::const_iterator i=simhits.begin(); i!= simhits.end(); ++i) {
      collectIonization(pixcharges, *i);
    }
  }

  //================================================================
  void ExtMonFNALHitMaker::collectIonization(PixelChargeCollection *pixcharges,
                                             const ExtMonFNALSimHit& hit)
  {
    const double sensorHalfThickness = extMon_->sensor().halfSize()[2];
    const double driftSpeed = siProps_.electronDriftMobility() * siProps_.electricField();


    // We split SimHit into a number of charge clusters, and drift them individually
    const CLHEP::Hep3Vector step = (1./(nclusters_-1)) * (hit.localEndPosition() - hit.localStartPosition());
    const double tstep = (1./(nclusters_-1)) * (hit.endTime() - hit.startTime());

    // charge in electrons
    const double meanClusterCharge = (1./nclusters_) * hit.ionizingEnergyDeposit() * siProps_.electronHolePairsPerEnergy();

    for(unsigned icluster=0; icluster<nclusters_; ++icluster) {
      const CLHEP::Hep3Vector pos = icluster*step + hit.localStartPosition();

      // The readout side is at z_local = -(sensor half thickness)
      // Floating point rounding can give negative values here, protect with max()
      const double driftDistance = std::max(0., pos.z() + sensorHalfThickness);

      const double driftTime = driftDistance/driftSpeed;

      const double diffusionSigma = std::sqrt(2.*siProps_.electronDiffusionConstant() * driftTime);

      // Readout positions
      const double x_ro = pos.x() + diffusionSigma * gaussian_.shoot();
      const double y_ro = pos.y() + diffusionSigma * gaussian_.shoot();

      // Charge fluctuations
      const double clusterCharge = meanClusterCharge +
        gaussian_.shoot() * std::sqrt(siProps_.fanoFactor() * meanClusterCharge);

      if(clusterCharge > 0) {
        // Record the contributions

        // While current flow starts at the ionization time, most of the pulse comes
        // from the pixel vicinity.  Roughly account for this by adding driftTime.
        const double time = hit.startTime() + icluster*tstep + driftTime;

        addCharge(pixcharges, hit.sensorId(), time, clusterCharge, x_ro, y_ro, hit.simParticle());
      }
    }
  }

  //================================================================
  void  ExtMonFNALHitMaker::addCharge(PixelChargeCollection *pixcharges,
                                      ExtMonFNALSensorId sid,
                                      double time,
                                      double charge,
                                      double x_ro,
                                      double y_ro,
                                      const art::Ptr<SimParticle>& particle)
  {
    ExtMonFNALPixelId pix = extMon_->sensor().findPixel(sid, x_ro, y_ro);
    if(pix != ExtMonFNALPixelId()) {
      (*pixcharges)[pix]
        .push(PixelTimedChargeDeposit(time, charge, particle));
    }
  }

  //================================================================
  void ExtMonFNALHitMaker::discriminate(ExtMonFNALRawHitCollection *outhits,
                                        ExtMonFNALHitTruthAssn *outtruth,
                                        art::ProductID hitsPID,
                                        const art::EDProductGetter *hitsGetter,
                                        PixelChargeCollection& pixcharges)
  {
    for(PixelChargeCollection::iterator i=pixcharges.begin(); i!=pixcharges.end(); ++i) {
      discriminate(outhits, outtruth, hitsPID, hitsGetter, i->first, i->second);
    }
  }

  //================================================================
  void ExtMonFNALHitMaker::discriminate(ExtMonFNALRawHitCollection *outhits,
                                        ExtMonFNALHitTruthAssn *outtruth,
                                        art::ProductID hitsPID,
                                        const art::EDProductGetter *hitsGetter,
                                        const ExtMonFNALPixelId& pix,
                                        PixelChargeHistory& ch)
  {
    PixelToTCircuit cap(discriminatorThreshold_, qCalib_, totCalib_, cond_->clockTick());

    double t = ch.top().time;

    while(!ch.empty()) {

      cap.wait(ch.top().time - t);
      t = ch.top().time;

      cap.addCharge(ch.top().charge);

      //----------------------------------------------------------------
      if(cap.high()) { // Found LE

        const int roStartTime = timeStamp(t);

        // add to the set of SimParticles
        typedef std::map<art::Ptr<SimParticle>, double> ChargeMap;
        ChargeMap parts;
        parts[ch.top().particle] += ch.top().charge;

        ch.pop();

        //----------------------------------------------------------------
        // Merge charges that go in the same hit

        while(!ch.empty() && (timeStamp(ch.top().time) <= timeStamp(t + cap.computeTrailingEdge()))) {

          cap.wait(ch.top().time - t);
          t = ch.top().time;

          cap.addCharge(ch.top().charge);
          parts[ch.top().particle] += ch.top().charge;

          ch.pop();
        }

        //----------------------------------------------------------------
        // Figure out ToT

        const int roEndTime = timeStamp(t + cap.computeTrailingEdge());

        // Limit dinamic range of readout ToT
        const int roToT = std::min(maxToT_, roEndTime - roStartTime);

        //----------------------------------------------------------------
        // Record the hit

        outhits->push_back(ExtMonFNALRawHit(pix, roStartTime, roToT));

        // Record hit truth
        for(ChargeMap::const_iterator t = parts.begin(); t != parts.end(); ++t) {

          outtruth->addSingle(t->first,

                              art::Ptr<ExtMonFNALRawHit>(hitsPID,
                                                         outhits->size()-1,
                                                         hitsGetter),

                              ExtMonFNALHitTruthBits(t->second)
                              );
        }

      }
      //----------------------------------------------------------------
      else {  // did not exceed the threshold, go on to the next charge
        ch.pop();
      }
    } // while(!empty)

  } // discriminate(pixel)

} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::ExtMonFNALHitMaker)
