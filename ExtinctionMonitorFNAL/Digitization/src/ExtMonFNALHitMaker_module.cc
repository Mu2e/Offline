// Pixel digitization: create ExtMonFNALRawHits and associated truth.
// Time stamps of created hits are in [0, numClockTicksPerDebuncherPeriod-1].
//
// Original author Andrei Gaponenko

#include <string>
#include <cmath>
#include <memory>
#include <queue>
#include <iostream>
#include <iomanip>

#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/Ptr.h"

#include "RecoDataProducts/inc/ExtMonFNALRawHit.hh"
#include "RecoDataProducts/inc/ExtMonFNALRawHitCollection.hh"
#include "MCDataProducts/inc/ExtMonFNALHitTruthAssn.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/ExtMonFNALSimHit.hh"
#include "MCDataProducts/inc/ExtMonFNALSimHitCollection.hh"

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALModule.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALModuleIdConverter.hh"
#include "ExtinctionMonitorFNAL/Digitization/inc/SiliconProperties.hh"
#include "ExtinctionMonitorFNAL/Digitization/inc/PixelToTCircuit.hh"
#include "ExtinctionMonitorFNAL/Digitization/inc/PixelNoise.hh"
#include "ExtinctionMonitorFNAL/Digitization/inc/PixelCharge.hh"
#include "ExtinctionMonitorFNAL/Digitization/inc/ProtonPulseShape.hh"

#include "GeometryService/inc/GeomHandle.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/ExtMonFNALConditions.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "SeedService/inc/SeedService.hh"


namespace mu2e {
  namespace ExtMonFNAL {

    struct VerilogHit {
      int pixelAddress;
      double twalk;
      double tot;
      VerilogHit(int p, double tw, double width) : pixelAddress(p), twalk(tw), tot(width) {}
    };

    struct VerilogHitAddrCmp {
      bool operator()(const VerilogHit& a, const VerilogHit& b) {
        return a.pixelAddress < b.pixelAddress;
      }
    };

    typedef std::vector<VerilogHit>  VerilogHitCollection;
    // BX => hits
    typedef std::map<int, VerilogHitCollection>  VerilogHitMap;

    //================================================================
    class ExtMonFNALHitMaker : public art::EDProducer {

    public:
      explicit ExtMonFNALHitMaker(fhicl::ParameterSet const& pset)
        : EDProducer{pset}
        , inputModuleLabel_(pset.get<std::string>("inputModuleLabel"))
        , inputInstanceName_(pset.get<std::string>("inputInstanceName"))
        , geomModuleLabel_(pset.get<std::string>("geomModuleLabel"))
        , geomInstanceName_(pset.get<std::string>("geomInstanceName", ""))
        , t0_(pset.get<double>("t0"))
        , nclusters_(pset.get<unsigned>("numClustersPerHit"))
        , maxToT_(pset.get<unsigned>("maxToT"))
        , discriminatorThreshold_(pset.get<double>("discriminatorThreshold"))
        , qCalib_(pset.get<double>("qCalib"))
        , totCalib_(pset.get<int>("totCalib"))

        , eng_(createEngine(art::ServiceHandle<SeedService>()->getSeed()))

        , gaussian_(eng_)

        , extMon_(0)
        , condExtMon_(0)
        , condAcc_(0)

        , applyProtonPulseShape_(pset.get<bool>("applyProtonPulseShape"))
        , protonPulse_(applyProtonPulseShape_ ?
                       new ProtonPulseShape(pset.get<fhicl::ParameterSet>("protonPulse"), eng_) :
                       0)

        , noise_(eng_,
                 &extMon_,/*Geometry not available at module ctr, store the address of the ptr */
                 &condExtMon_, /*similar for Conditions*/
                 pset.get<double>("pixelNoisePerBC"))

        , cutClockEnabled_(false)
        , cutClockMin_()
        , cutClockMax_()

        , chipSimInputsMode_(false)
        , chipSimFile_(nullptr)
        , chipSimChipId_()
        , chipSimProtonPulseNumber_()
      {
        produces<ExtMonFNALRawHitCollection>();
        produces<ExtMonFNALHitTruthAssn>();

        if(nclusters_ < 2) {
          throw cet::exception("CONFIG")
            <<"ExtMonFNALHitMaker: the numClustersPerHit FHICL parameter should be greater than 1\n";
        }

        std::cout<<"ExtMonFNALHitMaker: discriminatorThreshold = "<<discriminatorThreshold_
                 <<", qCalib = "<<qCalib_<<", totCalib = "<<totCalib_
                 <<std::endl;

        std::cout<<"ExtMonFNALHitMaker: t0 = "<<t0_<<", nclusters = "<<nclusters_<<std::endl;

        //----------------------------------------------------------------
        fhicl::ParameterSet cutClockPset;
        cutClockEnabled_ = pset.get_if_present("cutClock", cutClockPset);
        if(cutClockEnabled_) {
          cutClockMin_ = cutClockPset.get<int>("min");
          cutClockMax_ = cutClockPset.get<int>("max");
          std::cout<<"ExtMonFNALHitMaker: cutClock requested. Hits with clocks in ["
                   <<cutClockMin_<<", "<<cutClockMax_
                   <<"] will not be written out."<<std::endl;
        }
        else {
          std::cout<<"ExtMonFNALHitMaker: cutClock disabled"<<std::endl;
        }

        //----------------------------------------------------------------
        fhicl::ParameterSet chipSimPset;
        chipSimInputsMode_ = pset.get_if_present("chipSimInputs", chipSimPset);
        if(chipSimInputsMode_) {
          const std::string chipFileName(chipSimPset.get<std::string>("fileName"));
          chipSimFile_.reset(new std::ofstream(chipFileName.c_str()));
          ExtMonFNALModuleIdConverter con(*extMon_);
          ExtMonFNALModuleId mod = con.moduleId(ExtMonFNALModuleDenseId(chipSimPset.get<int>("module")));

          chipSimChipId_ = ExtMonFNALChipId(mod,
                                            chipSimPset.get<unsigned>("chipCol"),
                                            chipSimPset.get<unsigned>("chipRow")
                                            );

          std::cout<<"ExtMonFNALHitMaker: writing out Verilog inputs to file "<<chipFileName<<std::endl;
        }
      }

      virtual void produce(art::Event& evt);
      virtual void beginRun(art::Run& run);

    private:
      std::string inputModuleLabel_;
      std::string inputInstanceName_;
      std::string geomModuleLabel_;
      std::string geomInstanceName_;

      // global time corresponding to the starting edge of time bin 0
      double t0_;

      // The number of charge clusters per one simhit
      unsigned nclusters_;

      // Dynamic range of readout ToT
      int maxToT_;

      double discriminatorThreshold_;
      double qCalib_;
      int    totCalib_;

      art::RandomNumberGenerator::base_engine_t& eng_;
      CLHEP::RandGaussQ gaussian_;

      // Non-owning pointers to the geometry and conditions objects. The
      // current Mu2e infrastructure does not allow the use of a Handle
      // as a class member.
      const ExtMon *extMon_;
      const ExtMonFNALConditions *condExtMon_;
      const AcceleratorParams *condAcc_;

      SiliconProperties siProps_;

      bool  applyProtonPulseShape_;
      std::unique_ptr<ProtonPulseShape> protonPulse_;

      PixelNoise noise_;

      std::vector<double> planeTOFCorrection_;

      void collectIonization(PixelChargeCollection *pixcharges,
                             const ExtMonFNALSimHitCollection& simhits);

      void collectIonization(PixelChargeCollection *pixcharges,
                             const ExtMonFNALSimHit& hit);

      void addCharge(PixelChargeCollection *pixcharges,
                     ExtMonFNALModuleId mid,
                     double time,
                     double charge,
                     double x_ro,
                     double y_ro,
                     const art::Ptr<SimParticle>& particle);

      void foldHitTimes(PixelChargeCollection *inout);
      void foldHitTimes(PixelChargeHistory *inout);

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

      double hitTime_ns(unsigned iplane, double time) const {
        return (time - t0_ - planeTOFCorrection_[iplane]);
      }

      int timeStamp(double time_ns) const {
        return time_ns/condExtMon_->clockTick();
      }

      int timeStamp(int iplane, double time) const {
        return timeStamp(hitTime_ns(iplane, time));
      }

      //----------------
      // The cut on hit clock could be done by an external module or a consumer.
      // It is here for convenience.

      bool cutClockEnabled_;
      int cutClockMin_;
      int cutClockMax_;
      bool cutClockPassed(int clock) {
        return ! ((cutClockMin_ <= clock)&&(clock <= cutClockMax_));
      }

      //----------------
      // Write out a text file to be used as an input for the Verilog simulation
      bool chipSimInputsMode_;
      std::unique_ptr<std::ostream> chipSimFile_;
      ExtMonFNALChipId chipSimChipId_; // write out info for this single chip
      int chipSimProtonPulseNumber_;
      VerilogHitMap vlhm_;

      void addVerilogHit(const ExtMonFNALPixelId& pix, double tstart, double tend);
      void writeOutVerilogHits();

      //----------------
    };

    //================================================================
    void ExtMonFNALHitMaker::beginRun(art::Run& run) {
      if(!geomModuleLabel_.empty()) {
        art::Handle<ExtMon> emf;
        run.getByLabel(geomModuleLabel_, geomInstanceName_, emf);
        extMon_ = &*emf;
      }
      else {
        GeomHandle<ExtMon> emf;
        extMon_ = &*emf;
      }

      { // Got geometry. Compute per-plane time of flight correction to T0.
        const double p = extMon_->spectrometerMagnet().nominalMomentum();

        GlobalConstantsHandle<ParticleDataTable> pdt;
        ParticleDataTable::maybe_ref protonInfo = pdt->particle(2212);
        const double m = protonInfo.ref().mass();
        const double pm2 = std::pow(p/m, 2);
        const double beta = sqrt(pm2/(1.+pm2));
        const double v = beta * CLHEP::c_light;

        planeTOFCorrection_.resize(extMon_->nplanes());
        for(unsigned i=0; i<extMon_->nplanes(); ++i) {
          const CLHEP::Hep3Vector pos = extMon_->planeCenterInExtMon(i);
          planeTOFCorrection_[i] = -pos.z()/v;
        }
      }

      ConditionsHandle<ExtMonFNALConditions> cond("ignored");
      condExtMon_ = &*cond;

      ConditionsHandle<AcceleratorParams> condAcc("ignored");
      condAcc_ = &*condAcc;

      const double moduleThickness = 2*(extMon_->module().sensorHalfSize()[2] + extMon_->module().chipHalfSize()[2]);

      const double electricField = std::abs(cond->biasVoltage())/moduleThickness;

      siProps_.setConditions(cond->temperature(), electricField);

      std::cout<<"ExtMonFNALHitMaker:"
               <<"  electronHolePairsPerEnergy = "<<siProps_.electronHolePairsPerEnergy()
               <<", electronDriftMobility = "<<siProps_.electronDriftMobility()
               <<", electronDiffusionConstant = "<<siProps_.electronDiffusionConstant()
               <<std::endl;

      std::cout<<"ExtMonFNALHitMaker:"
               <<"  numClockTicksPerDebuncherPeriod = "<<cond->numClockTicksPerDebuncherPeriod()
               <<", clockTick = "<<cond->clockTick()
               <<std::endl;
    }

    //================================================================
    void ExtMonFNALHitMaker::produce(art::Event& event) {
      std::unique_ptr<ExtMonFNALRawHitCollection> outHits(new ExtMonFNALRawHitCollection());
      std::unique_ptr<ExtMonFNALHitTruthAssn> outTruth(new ExtMonFNALHitTruthAssn());

      art::Handle<ExtMonFNALSimHitCollection> ih;
      event.getByLabel(inputModuleLabel_, inputInstanceName_, ih);
      const ExtMonFNALSimHitCollection& simhits(*ih);

      chipSimProtonPulseNumber_ = event.event();

      // N.B.: to reduce memory footprint we can digitize one readout chip at a time
      // instead of doing the whole detector in one go.  Would be an extra loop here.
      PixelChargeCollection pixcharges;
      collectIonization(&pixcharges, simhits);

      if(applyProtonPulseShape_) {
        protonPulse_->apply(&pixcharges, event);
      }

      // Brings all times onto a microbunch + margins on both sides.
      // Hits near microbunch boundaries are duplicated.
      foldHitTimes(&pixcharges);

      const art::ProductID hitsPID = event.getProductID<ExtMonFNALRawHitCollection>();
      const art::EDProductGetter *hitsGetter = event.productGetter(hitsPID);
      discriminate(&*outHits, &*outTruth, hitsPID, hitsGetter, pixcharges);

      if(chipSimInputsMode_) {
        writeOutVerilogHits();
        vlhm_.clear();
      }

      noise_.add(&*outHits);

      event.put(std::move(outHits));
      event.put(std::move(outTruth));
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
      const double moduleHalfThickness = (extMon_->module().sensorHalfSize()[2] + extMon_->module().chipHalfSize()[2]);
      const double driftSpeed = siProps_.electronDriftMobility() * siProps_.electricField();


      // We split SimHit into a number of charge clusters, and drift them individually
      const CLHEP::Hep3Vector step = (1./(nclusters_-1)) * (hit.localEndPosition() - hit.localStartPosition());
      const double tstep = (1./(nclusters_-1)) * (hit.endTime() - hit.startTime());

      // Average deposited charge in electrons
      const double meanTotalCharge = hit.ionizingEnergyDeposit() * siProps_.electronHolePairsPerEnergy();

      // Charge fluctuations
      const double totalCharge = meanTotalCharge +
        gaussian_.fire() * std::sqrt(siProps_.fanoFactor() * meanTotalCharge);

      const double clusterCharge = totalCharge/nclusters_;

      if(clusterCharge > 0) { // there are occasional energy deposits corresponding to < 1 pair, which can fluctuate to <0
        for(unsigned icluster=0; icluster<nclusters_; ++icluster) {
          const CLHEP::Hep3Vector pos = icluster*step + hit.localStartPosition();

          // The readout side is at z_local = -(module half thickness)
          // Floating point rounding can give negative values here, protect with max()
          const double driftDistance = std::max(0., pos.z() + moduleHalfThickness);

          const double driftTime = driftDistance/driftSpeed;

          const double diffusionSigma = std::sqrt(2.*siProps_.electronDiffusionConstant() * driftTime);

          // Readout positions
          const double x_ro = pos.x() + diffusionSigma * gaussian_.fire();
          const double y_ro = pos.y() + diffusionSigma * gaussian_.fire();

          // While current flow starts at the ionization time, most of the pulse comes
          // from the pixel vicinity.  Roughly account for this by adding driftTime.
          const double time = hit.startTime() + icluster*tstep + driftTime;

          addCharge(pixcharges, hit.moduleId(), time, clusterCharge, x_ro, y_ro, hit.simParticle());
        }
      }
    }

    //================================================================
    void  ExtMonFNALHitMaker::addCharge(PixelChargeCollection *pixcharges,
                                        ExtMonFNALModuleId mid,
                                        double time,
                                        double charge,
                                        double x_ro,
                                        double y_ro,
                                        const art::Ptr<SimParticle>& particle)
    {
      ExtMonFNALPixelId pix = extMon_->module().findPixel(mid, x_ro, y_ro);
      if(pix != ExtMonFNALPixelId()) {
        (*pixcharges)[pix]
          .push(PixelTimedChargeDeposit(time, charge, particle));
      }
    }

    //================================================================
    void ExtMonFNALHitMaker::foldHitTimes(PixelChargeCollection* inout) {
      for(PixelChargeCollection::iterator i=inout->begin(); i!=inout->end(); ++i) {
        foldHitTimes(&i->second);
      }
    }

    //================================================================
    void ExtMonFNALHitMaker::foldHitTimes(PixelChargeHistory* inout) {
      // Here we bring hits from (-infty, +infty) to
      // (-margin, deBuncherPeriod + margin)
      //
      // Margins on both sides make sure hits near both of the
      // microbunch boundaries are modeled correctly.  Hits in +-
      // margin within a boundary are duplicated, then one of them is
      // cut off after discrimination in the final "time folding"
      // step, which is truncation of discrete hit times.  (Hits
      // ending up at t=-delta<0 do not produce output hits, but may
      // eat up other hits with t>0, whose effect is instead in
      // modifying ToT of the twin hit t=deBuncherPeriod-delta.)
      //
      // Larger margins are safe (for correctnes), but inefficient.
      // maxToT + max time of flight correction should be enough.


      PixelChargeHistory& in(*inout);
      PixelChargeHistory out;

      const double period = condAcc_->deBuncherPeriod;
      const double margin = maxToT_ * condExtMon_->clockTick();

      while(!in.empty()) {
        PixelTimedChargeDeposit dep = in.top();
        // A hit on [0, period]
        dep.time = remainder(dep.time - period/2, period) + period/2;
        out.push(dep);

        // duplicate hits near the boundaries
        if(dep.time < margin) {
          dep.time = dep.time + period;
          out.push(dep);
        }
        else if(period - dep.time < margin) {
          dep.time = dep.time - period;
          out.push(dep);
        }

        in.pop();
      }

      *inout = out;
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
      PixelToTCircuit cap(discriminatorThreshold_, qCalib_, totCalib_, condExtMon_->clockTick());
      const unsigned iplane = pix.chip().module().plane();

      double t = ch.top().time;

      while(!ch.empty()) {

        cap.wait(ch.top().time - t);
        t = ch.top().time;

        cap.addCharge(ch.top().charge);

        //----------------------------------------------------------------
        if(cap.high()) { // Found LE

          const double hitStart_ns = hitTime_ns(iplane, t);
          const int roStartTime = timeStamp(hitStart_ns);

          // add to the set of SimParticles
          typedef std::map<art::Ptr<SimParticle>, double> ChargeMap;
          ChargeMap parts;
          parts[ch.top().particle] += ch.top().charge;

          ch.pop();

          //----------------------------------------------------------------
          // Merge charges that go in the same hit

          while(!ch.empty() && (timeStamp(iplane, ch.top().time) <= timeStamp(iplane, t + cap.computeTrailingEdge()))) {

            cap.wait(ch.top().time - t);
            t = ch.top().time;

            cap.addCharge(ch.top().charge);
            parts[ch.top().particle] += ch.top().charge;

            ch.pop();
          }

          //----------------------------------------------------------------
          // Figure out ToT

          // NB: for chip sim we could have stopped the above loop
          // earlier.  However let's keep it as is to reduce effects
          // of the artificial "sawtooth" charge collection shape.
          const double hitEnd_ns = hitTime_ns(iplane, t + cap.computeTrailingEdge());
          const int roEndTime = timeStamp(hitEnd_ns);

          // Limit dinamic range of readout ToT
          const int roToT = std::min(maxToT_, roEndTime - roStartTime);

          //----------------------------------------------------------------
          // Record the hit
          //
          // For the proper time folding we should truncate here to
          // one microbunch, which compensates for hit duplication in
          // analogue time folding.

          if((0 <= roStartTime) && (roStartTime < condExtMon_->numClockTicksPerDebuncherPeriod())) {

            if(chipSimInputsMode_) {
              addVerilogHit(pix, hitStart_ns, hitEnd_ns);
            }

            if(!cutClockEnabled_ || cutClockPassed(roStartTime)) {
              outhits->push_back(ExtMonFNALRawHit(pix, roStartTime, roToT));

              // Record hit truth
              for(ChargeMap::const_iterator t = parts.begin(); t != parts.end(); ++t) {

                outtruth->addSingle(t->first,

                                    art::Ptr<ExtMonFNALRawHit>(hitsPID,
                                                               outhits->size()-1,
                                                               hitsGetter),

                                    ExtMonFNALHitTruthBits(t->second)
                                    );
              } // for()

            } // cutClock
          } // truncation to debuncher period

        }
        //----------------------------------------------------------------
        else {  // did not exceed the threshold, go on to the next charge
          ch.pop();
        }
      } // while(!empty)

    } // discriminate(pixel)


    //================================================================
    void ExtMonFNALHitMaker::addVerilogHit(const ExtMonFNALPixelId& pix, double tstart, double tend) {
      if(pix.chip() == chipSimChipId_) {
        const int bx = timeStamp(tstart) +
          condExtMon_->numClockTicksPerDebuncherPeriod() * (chipSimProtonPulseNumber_ - 1);

        const double twalk = tstart - timeStamp(tstart) * condExtMon_->clockTick();

        const double tot = tend - tstart;

        //----------------------------------------------------------------
        // Date: Thu, 21 Feb 2013 00:02:04
        // From: Dario Gnani <dgnani@lbl.gov>
        // To: Andrei Gaponenko <gandr@fnal.gov>
        // Cc: maurice <mgarcia-sciveres@lbl.gov>, David Brown <Dave_Brown@lbl.gov>
        // Subject: Re: FE-I4 sim for Mu2e
        //
        // [...]
        //
        // the pixel address encoding for (col,row) with col in [1:80] and row in [1:336] is
        // pix_addr=336*(col-1)+(row-1)
        // Physically -looking at the ASIC with the pads at the bottom, column
        // 1 is on the left and row 1 at the top
        //
        // AG: in Mu2e Offline we use 0 based pixel col and row, so:

        const int pixelAddress = 336 * pix.col() + pix.row();

        vlhm_[bx].push_back(VerilogHit(pixelAddress, twalk, tot));
      }
    }

    //================================================================
    void ExtMonFNALHitMaker::writeOutVerilogHits() {

      // The OVM input script wants to have a line for each clock tick
      // Therefore we do a fixed-size loop here and write out a lot of emtpy events.
      for(int localClock = 0; localClock < condExtMon_->numClockTicksPerDebuncherPeriod(); ++localClock) {

        const int ibx =  localClock + condExtMon_->numClockTicksPerDebuncherPeriod() * (chipSimProtonPulseNumber_ - 1);

        *chipSimFile_ <<"BX "<<ibx<<std::endl;

        VerilogHitMap::const_iterator pd = vlhm_.find(ibx);
        if(pd != vlhm_.end()) {

          //NB: here we make a copy to apply non-const sort
          VerilogHitCollection coll(pd->second);
          std::sort(coll.begin(), coll.end(), VerilogHitAddrCmp());
          for(VerilogHitCollection::const_iterator i = coll.begin(); i != coll.end(); ++i) {
            *chipSimFile_<<i->pixelAddress
                         <<"\t"
                         <<std::fixed<<std::setprecision(0)
                         <<i->twalk
                         <<"\t"<<i->tot
                         <<std::endl;
          }
        }
      }
    }

    //================================================================

  } // end namespace ExtMonFNAL
} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::ExtMonFNAL::ExtMonFNALHitMaker)
