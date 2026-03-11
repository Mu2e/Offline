#include "Offline/CaloReco/inc/CaloTemplateWFProcessor.hh"
#include "Offline/CaloReco/inc/CaloTemplateWFUtil.hh"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

#include "TFile.h"
#include "TH2.h"

#include <algorithm>
#include <string>
#include <iostream>
#include <vector>
#include <numeric>

namespace mu2e {

   CaloTemplateWFProcessor::CaloTemplateWFProcessor(const Config& config) :
      CaloWaveformProcessor(),
      windowPeak_      (config.windowPeak()),
      minPeakAmplitude_(config.minPeakAmplitude()),
      numNoiseBins_    (config.numNoiseBins()),
      minDTPeaks_      (config.minDTPeaks()),
      psdThreshold_    (config.psdThreshold()),
      chiThreshold_    (config.chiThreshold()),
      refitLeadingEdge_(config.refitLeadingEdge()),
      diagLevel_       (config.diagLevel()),
      fmutil_          (config.pulseFileName(),config.pulseHistName(),minPeakAmplitude_,
                        config.digiSampling(),minDTPeaks_,config.fitPrintLevel()),
      chi2_            (999.),
      ndf_             (-1),
      resAmp_          (),
      resAmpErr_       (),
      resTime_         (),
      resTimeErr_      ()
   {
       if (diagLevel_ > 1) initHistos();
       if (windowPeak_ < 1) windowPeak_=1;
   }



   void CaloTemplateWFProcessor::initialize()
   {
       fmutil_.initialize();
   }


   void CaloTemplateWFProcessor::extract(const std::vector<double>& xInput, const std::vector<double>& yInput)
   {
       if (diagLevel_>2) std::cout<<"CaloTemplateWFProcessor start"<<std::endl;

       reset();
       fmutil_.setXYVector(xInput, yInput);


       //try first strategy
       setPrimaryPeakPar1(xInput, yInput);
       if (fmutil_.nPeaks()==0) return;
       setSecondaryPeakPar(xInput, yInput);

       fmutil_.fit();


       //if that fails, try a more complicated model
       chi2_ = fmutil_.chi2();
       ndf_  = yInput.size() - fmutil_.par().size();
       if (ndf_ > 0 && chi2_/float(ndf_)>chiThreshold_)
       {
            setPrimaryPeakPar2(xInput, yInput);
            fmutil_.fit();

            if (fmutil_.maxAmplitude() > 3*minPeakAmplitude_)
            {
                auto nFirstPeaks = fmutil_.nPeaks();
                setSecondaryPeakPar(xInput, yInput);
                if (fmutil_.nPeaks() > nFirstPeaks) fmutil_.fit();
            }
       }


       if (refitLeadingEdge_) fmutil_.refitEdge();

       for (unsigned i=0;i<fmutil_.nPeaks();++i)
       {
           unsigned ip = fmutil_.peakIdx(i);
           if (fmutil_.par()[ip+1] > xInput.back()) continue;
           resAmp_.push_back(fmutil_.par()[ip]);
           //modified estimate for the uncertainty, see note in CaloTemplateWFUtil
           resAmpErr_.push_back( sqrt(fmutil_.parErr()[ip]*fmutil_.parErr()[ip]+fmutil_.par()[ip]) );
           resTime_.push_back(fmutil_.fromPeakToT0(fmutil_.par()[ip+1]));
           resTimeErr_.push_back(fmutil_.parErr()[ip+1]);
       }

       chi2_ = fmutil_.chi2();
       ndf_  = yInput.size() - fmutil_.par().size();

       if (diagLevel_ > 1)
       {
          double chindf = (ndf_ > 0) ? chi2_/float(ndf_) : 999;
          _hNpeak->Fill(resAmp_.size());
          _hChi2->Fill(chindf);
          _hchi2Peak->Fill(chindf,resAmp_.size());
          for (auto val : resAmp_)     _hchi2Amp->Fill(chindf,val);
          for (auto val : resAmp_)     _hEner->Fill(val);
          for (auto val : resAmpErr_)  _hEnerErr->Fill(val);
          for (auto val : resTime_)    _hTime->Fill(val);
          for (auto val : resTimeErr_) _hTimeErr->Fill(val);
       }

       if (diagLevel_>2) std::cout<<"Fit done"<<std::endl;
   }



   //----------------------------------------------------------------------------------------------------------------------------------
   void CaloTemplateWFProcessor::reset() {resAmp_.clear(); resAmpErr_.clear(); resTime_.clear(); resTimeErr_.clear();ndf_ = 0;chi2_ = 999; fmutil_.reset();}

   //---------------------------------------------------------------------------------------------------------------------------------------------
   void CaloTemplateWFProcessor::setPrimaryPeakPar1(const std::vector<double>& xvec, const std::vector<double>& yvecOrig)
   {
        if (windowPeak_ > xvec.size()) return;
        std::vector<double> parInit{},yvec(yvecOrig);

        //estimate the noise level with the first few bins
        if (yvec.size() <= numNoiseBins_) return;
        float noise= std::accumulate(yvec.begin(),yvec.begin()+numNoiseBins_,0)/float(numNoiseBins_);
        parInit.push_back(noise);
        for (auto& val : yvec) val -= noise;

        if (diagLevel_>1) std::cout<<"[CaloTemplateWFProcessor] Noise level "<<noise<<std::endl;

        std::vector<double> ywork(yvec);
        for (auto i=windowPeak_;i+windowPeak_<xvec.size();++i)
        {
            if (std::max_element(ywork.begin()+i-windowPeak_,ywork.begin()+i+windowPeak_+1) != ywork.begin()+i) continue;
            if (ywork[i-1] < minPeakAmplitude_ || ywork[i] < minPeakAmplitude_ || ywork[i+1] < minPeakAmplitude_) continue;

            double xmaxEstimate = estimatePeakTime(xvec,ywork,i);
            double ymaxEstimate = fmutil_.peakNorm(xvec, ywork, xmaxEstimate , i-windowPeak_, i+windowPeak_);

            parInit.push_back(ymaxEstimate);
            parInit.push_back(xmaxEstimate);
            fmutil_.setPar(parInit);

            if (diagLevel_>1) std::cout<<"[CaloTemplateWFProcessor] Found primary peak "<<xmaxEstimate<<std::endl;
        }

        if (fmutil_.nPeaks()==0)
        {
            unsigned ipeak(0);
            while (ipeak+1 < ywork.size()) {if (ywork[ipeak] >= minPeakAmplitude_ && ywork[ipeak+1] >= 0.8*ywork[ipeak]) break; ++ipeak;}
            if (ipeak+1 < ywork.size())
            {
                parInit.push_back(ywork[ipeak]);
                parInit.push_back(xvec[ipeak]);
                fmutil_.setPar(parInit);
                if (diagLevel_>1) std::cout<<"[CaloTemplateWFProcessor] Found primary peak alt "<<xvec[ipeak] <<std::endl;
            }
        }

    }


   //---------------------------------------------------------------------------------------------------------------------------------------------
   void CaloTemplateWFProcessor::setPrimaryPeakPar2(const std::vector<double>& xvec, const std::vector<double>& yvecOrig)
   {
        if (windowPeak_ > xvec.size()) return;
        std::vector<double> parInit{},yvec(yvecOrig);

        //estimate the noise level with the first few bins
        float noise= std::accumulate(yvec.begin(),yvec.begin()+numNoiseBins_,0)/float(numNoiseBins_);
        parInit.push_back(noise);
        for (auto& val : yvec) val -= noise;

        std::vector<double> ywork(yvec);
        for (auto i=windowPeak_;i<int(xvec.size())-windowPeak_;++i)
        {
            if (std::max_element(ywork.begin()+i-windowPeak_,ywork.begin()+i+windowPeak_+1) != ywork.begin()+i) continue;
            if (ywork[i-1] < minPeakAmplitude_ || ywork[i] < minPeakAmplitude_ || ywork[i+1] < minPeakAmplitude_) continue;

            double xmaxEstimate = estimatePeakTime(xvec,ywork,i);
            double ymaxEstimate = fmutil_.peakNorm(xvec, ywork, xmaxEstimate , i-windowPeak_, i+windowPeak_);

            parInit.push_back(ymaxEstimate);
            parInit.push_back(xmaxEstimate);

            if (diagLevel_ > 2) std::cout<<"[CaloTemplateWFProcessor] Found primary peak bis "<<xmaxEstimate<<std::endl;

            fmutil_.setPar(parInit);
            for (unsigned j=0;j<xvec.size();++j) ywork[j] -= fmutil_.eval_logn(xvec[j],parInit.size()-fmutil_.nParFcn());

            if (ymaxEstimate > 3*minPeakAmplitude_) findRisingPeak(i, parInit, xvec, yvec, ywork);
        }

        if (fmutil_.nPeaks()==0)
        {
            unsigned ipeak(0);
            while (ipeak < ywork.size()) {if (ywork[ipeak] >= minPeakAmplitude_) break; ++ipeak;}
            if (ipeak < ywork.size())
            {
                parInit.push_back(ywork[ipeak]);
                parInit.push_back(xvec[ipeak]);
                fmutil_.setPar(parInit);
            }
        }
    }

    //--------------------------------------------------------------------------------------------------
    void CaloTemplateWFProcessor::findRisingPeak(int ipeak, std::vector<double>& parInit, const std::vector<double>& xvec,
                                                 const std::vector<double>& yvec, std::vector<double>& ywork)
    {
        //find where the peak "starts" - check it is far away from previous peak
        unsigned istartPeak(ipeak - windowPeak_);
        while (istartPeak>windowPeak_ && yvec[istartPeak] > 0.5*yvec[ipeak]) --istartPeak;
        while (istartPeak>windowPeak_ && yvec[istartPeak] > minPeakAmplitude_ && yvec[istartPeak] > yvec[istartPeak-1]) --istartPeak;

        //find max of the residuals in the rising edge of peak
        auto peakMax     = std::max_element(ywork.begin()+istartPeak,ywork.begin()+ipeak-2);
        unsigned peakIdx = std::distance(std::begin(ywork),peakMax);
        double minAmp    =  ywork[ipeak]/3.0;

        if (ywork[peakIdx] > minAmp)
        {
            double tl(xvec[peakIdx]),th(parInit.back()),dt(0.2*(xvec[1]-xvec[0]));
            while (th - tl > dt)
            {
               if (fmutil_.sumSquare(xvec,ywork, tl, peakIdx-windowPeak_,peakIdx) >
                   fmutil_.sumSquare(xvec,ywork,th,peakIdx-windowPeak_,peakIdx)) tl = 0.5*(tl+th);
               else                                                              th = 0.5*(tl+th);
            }

            double peakTime  = tl;
            double Amplitude = fmutil_.peakNorm(xvec,yvec, peakTime, peakIdx-windowPeak_,peakIdx);

            if (Amplitude > minAmp)
            {
                for (unsigned j=0;j<xvec.size();++j) ywork[j] += fmutil_.eval_logn(xvec[j],parInit.size()-fmutil_.nParFcn());

                //replace old peak by new rising peak
                parInit[parInit.size()-2] = Amplitude;
                parInit[parInit.size()-1] = peakTime;

                fmutil_.setPar(parInit);
                for (unsigned j=0;j<xvec.size();++j) ywork[j] -= fmutil_.eval_logn(xvec[j],parInit.size()-fmutil_.nParFcn());
                if (diagLevel_ > 2) std::cout<<"Found rising peak "<<Amplitude<<" "<<peakTime<<std::endl;
            }
        }
    }



    //--------------------------------------------------------------------------------------------------
    void CaloTemplateWFProcessor::setSecondaryPeakPar(const std::vector<double>& xvec, const std::vector<double>& yvec)
    {
        if (windowPeak_ > xvec.size()) return;
        if (xvec.size() != yvec.size()) throw cet::exception("CATEGORY")<<"CaloTemplateWFProcessor::setSecondaryPeakPar  xvec and yvec must have the same size";

        std::vector<double> ywork(yvec);
        for (unsigned j=0;j<xvec.size();++j) ywork[j] -= fmutil_.eval_fcn(xvec[j]);

        std::vector<double> parInit(fmutil_.par());
        for (auto i=windowPeak_;i<ywork.size()-windowPeak_;++i)
        {
             if (std::max_element(ywork.begin()+i-windowPeak_,ywork.begin()+i+windowPeak_+1) != ywork.begin()+i) continue;
             if (ywork[i-1] < minPeakAmplitude_ || ywork[i] < minPeakAmplitude_ || ywork[i+1] < minPeakAmplitude_) continue;
             if (checkPeakDist(xvec[i])) continue;
             if (yvec[i]>0 && ywork[i]/yvec[i] < psdThreshold_) continue;

             double xmaxEstimate = estimatePeakTime(xvec,ywork,i);
             //double xmaxEstimate = xvec[i];

             parInit.push_back(ywork[i]);
             parInit.push_back(xmaxEstimate);

             fmutil_.setPar(parInit);
             for (unsigned j=0;j<xvec.size();++j) ywork[j] -= fmutil_.eval_logn(xvec[j],parInit.size()-fmutil_.nParFcn());

             if (diagLevel_>2)  std::cout<<"Secondary peak "<<xvec[i]<<" "<<xmaxEstimate<<std::endl;
        }
    }


   //------------------------------------------------------------
   double CaloTemplateWFProcessor::estimatePeakTime(const std::vector<double>& xvec, const std::vector<double>& ywork, int ic)
   {
       double y1      = fmutil_.sumSquare(xvec, ywork, xvec[ic]  , ic-windowPeak_, ic+windowPeak_);
       double y2      = fmutil_.sumSquare(xvec, ywork, xvec[ic]-1, ic-windowPeak_, ic+windowPeak_);
       double y3      = fmutil_.sumSquare(xvec, ywork, xvec[ic]-2, ic-windowPeak_, ic+windowPeak_);
       double a       = (y1-2*y2+y3)/2;
       double b       = (y1-y2-2*a*xvec[ic]+a);
       double newTime = (std::abs(a)> 1e-6) ? -b/2.0/a : 0;

       if (std::abs(newTime-xvec[ic])<15.0) return newTime;
       return  xvec[ic];
   }

   //------------------------------------------------------------------------------------------
   bool CaloTemplateWFProcessor::checkPeakDist(double x0)
   {
       for (unsigned i=0;i<fmutil_.nPeaks();++i)
       {
           unsigned ip = fmutil_.peakIdx(i);
           if (std::abs(fmutil_.par()[ip+1]-x0) < minDTPeaks_) return true;
       }
       return false;
   }

   //---------------------------------------------------------------------------------------------------------------------------------------
   void CaloTemplateWFProcessor::plot(const std::string& name) const {fmutil_.plotFit(name);}

   //---------------------------------------------------------------------------------------------------------------------------------------
   void CaloTemplateWFProcessor::dump(const std::string& name, const std::vector<double>& val) const
   {
      std::cout<<name<<std::endl;
      for (auto h : val) std::cout<<h<<" ";
      std::cout<<std::endl;
   }

   //---------------------------------------------------------------------------------------------------------------------------------------
   void CaloTemplateWFProcessor::initHistos()
   {
       art::ServiceHandle<art::TFileService> tfs;
       art::TFileDirectory tfdir = tfs->mkdir("CaloTemplateWFDiag");
       _hTime     = tfdir.make<TH1F>("hTime",    "time;t (ns); Entries",                  100, 0., 2000);
       _hTimeErr  = tfdir.make<TH1F>("hTimeErr", "time error;dt (ns); Entries",           100, 0.,   10);
       _hEner     = tfdir.make<TH1F>("hAmp",     "Amplitude;A;Entries",                   200, 0.,  200);
       _hEnerErr  = tfdir.make<TH1F>("hAmpErr",  "Amplitude error;A_{err};Entries",       100, 0.,  100);
       _hNpeak    = tfdir.make<TH1F>("hNpeak",   "Number of peak fitted;N_{peak};Entries", 10, 0.,   10);
       _hChi2     = tfdir.make<TH1F>("hChi2",    "Chi2/ndf;#chi^{2}/ndf;Entries",         100, 0.,   20);
       _hchi2Amp  = tfdir.make<TH2F>("hchi2Amp", "Amp vs chi2/ndf",    50, 0.,   20, 100, 0, 2000);
       _hchi2Peak = tfdir.make<TH2F>("hchi2Peak","Npeak vs chi2/ndf",  50, 0.,   20, 10, 0, 10);
   }

}


/*
   //alternative peak finding algorithm, performs slightly below the others but kept for the record
   void CaloTemplateWFProcessor::setPrimaryPeakPar1(const std::vector<double>& xvec, const std::vector<double>& yvec)
   {

        if (2*windowPeak_ >= xvec.size()) return;
        std::vector<double> parInit,ywork(yvec);

        //estimate the noise level with the first few bins
        float noise= std::accumulate(yvec.begin(),yvec.begin()+numNoiseBins_,0)/float(numNoiseBins_);
        for (auto& val : ywork) val -= noise;
        parInit.push_back(noise);

        unsigned maxWF = *std::max_element(ywork.begin()+windowPeak_,ywork.end()-windowPeak_);

        unsigned nTry(0);
        while (nTry <20)
        {
            std::vector<double> ywork2(ywork);
            for (unsigned i=windowPeak_; i<ywork2.size()-windowPeak_;++i) ywork2[i] = (ywork[i-1]+ywork[i]+ywork[i+1])/3.0;
            auto maxEl = std::max_element(ywork2.begin()+windowPeak_,ywork2.end()-windowPeak_);
            int i = std::distance(ywork2.begin(),maxEl);
            if (ywork2[i] < minPeakAmplitude_ || ywork2[i-1] < minPeakAmplitude_ || ywork2[i-1] < minPeakAmplitude_) break;

            double xmaxEstimate = estimatePeakTime(xvec, ywork, i);
            double ymaxEstimate = fmutil_.peakNorm(xvec, ywork, xmaxEstimate , i-windowPeak_, i+windowPeak_);

            if (diagLevel_ > 2) std::cout<<"Found peak "<<xmaxEstimate<<" "<<ymaxEstimate<<" "<<xvec[i]<<std::endl;

            if (ymaxEstimate < maxWF/3.0 && std::max_element(yvec.begin()+i-windowPeak_,yvec.begin()+i+windowPeak_) != yvec.begin()+i)
            {
                for (unsigned j=i-windowPeak_;j<i+windowPeak_;++j) ywork[j] =0;
                ++nTry;
                continue;
            }

            std::vector<double> tempPar{noise};
            for (unsigned i=0;i<fmutil_.nPeaks();++i)
            {
                unsigned ip = fmutil_.peakIdx(i);
                if (parInit[ip+1] > xmaxEstimate && parInit[ip+1] < xmaxEstimate+50.0 && fmutil_.peakToFunc(ip, xmaxEstimate, ymaxEstimate) > 0.1)
                {
                    for (unsigned j=0;j<xvec.size();++j) ywork[j] += fmutil_.eval_logn(xvec[j],ip);
                }
                else
                {
                   tempPar.push_back(parInit[ip]);
                   tempPar.push_back(parInit[ip+1]);
                }
            }

            parInit = tempPar;
            parInit.push_back(ymaxEstimate);
            parInit.push_back(xmaxEstimate);

            fmutil_.setPar(parInit);
            for (unsigned j=0;j<xvec.size();++j) ywork[j] -= fmutil_.eval_logn(xvec[j],parInit.size()-fmutil_.nParFcn());

            ++nTry;
        }

        if (fmutil_.nPeaks()==0)
        {
            ywork = yvec;

            unsigned ipeak(0);
            while (ipeak < ywork.size()) {if (ywork[ipeak] >= minPeakAmplitude_) break; ++ipeak;}
            if (ipeak < ywork.size())
            {
                parInit.push_back(ywork[ipeak]);
                parInit.push_back(xvec[ipeak]);
                fmutil_.setPar(parInit);
            }
        }

        if (diagLevel_ > 2)
        {
            for (unsigned i=0;i<fmutil_.nPeaks();++i)
            {
               unsigned ip = fmutil_.peakIdx(i);
               std::cout<<"found peak at t="<<fmutil_.par()[ip+1]<<std::endl;
            }
        }
     }
*/
