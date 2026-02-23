#include "cetlib_except/exception.h"
#include "Offline/Mu2eUtilities/inc/CaloPulseShape.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"

#include "TFile.h"
#include "TH2F.h"

#include <memory>
#include <vector>
#include <iostream>


namespace mu2e {

   CaloPulseShape::CaloPulseShape(const std::string& fileName, const std::string& histName, double digiSampling) :
      fileName_(fileName),
      histName_(histName),
      nSteps_(100),
      digiStep_(digiSampling/double(nSteps_)),
      nBinShape_(0),
      pulseVec_(),
      deltaT_(0.),
      digitizedPulse_()
   {}

   //----------------------------------------------------------------------------------------------------------------------
   void CaloPulseShape::buildShapes()
   {

       pulseVec_.clear();

       ConfigFileLookupPolicy resolveFullPath;
       std::string fullFileName = resolveFullPath(fileName_);

       std::unique_ptr<TH1F> pshape(nullptr);
       TFile pulseFile(fullFileName.c_str());
       if (pulseFile.IsOpen()) pshape.reset((TH1F*) pulseFile.Get(histName_.c_str()));
       if (!pshape) throw cet::exception("CATEGORY")<<"CaloPulseShape:: Hitsogram "<<histName_.c_str()
                                                    <<" from file "<<fileName_.c_str()<<" does not exist";
       pshape->SetDirectory(0);
       pulseFile.Close();


       // Adjust binning to match digitizer sampling period, shift to zero and normalize
       int nbins = int((pshape->GetXaxis()->GetXmax()-pshape->GetXaxis()->GetXmin())/digiStep_);
       TH1F pulseShape("ps","ps", nbins, 0.0, pshape->GetXaxis()->GetXmax()-pshape->GetXaxis()->GetXmin());
       for (int i=1;i<=nbins;++i) pulseShape.SetBinContent(i,pshape->Interpolate(pulseShape.GetBinCenter(i)));
       pulseShape.Scale(1.0/pulseShape.GetMaximum(),"nosw2");

       // Cache histogram content into vector and shift waveform (see note),
       // calculate the number of bins for the digitized waveform
       for (int j=1;j<=(nbins+nSteps_);++j)pulseVec_.push_back((j>nSteps_) ? pulseShape.GetBinContent(j-nSteps_) : 0.0);
       nBinShape_      = int(nbins/nSteps_);
       digitizedPulse_ = std::vector<double>(nBinShape_,0);

       deltaT_ = 0.0;
       // find difference between peak time and t0 for digitized waveform.
       for (int i=1;i<nBinShape_;++i) {if (pulseVec_[(i+1)*nSteps_] < pulseVec_[i*nSteps_]) break; deltaT_ +=nSteps_*digiStep_;}

   }

   //----------------------------------------------------------------------------
   // forward shift in waveform = backward shift in time origin
   const std::vector<double>& CaloPulseShape::digitizedPulse(double hitTime) const
   {
       int shiftBin = nSteps_ - int(hitTime/digiStep_)%nSteps_;
       for (int i=0;i<nBinShape_;++i) digitizedPulse_[i] = pulseVec_[shiftBin+i*nSteps_];
       return digitizedPulse_;
   }

   //----------------------------------------------------------------------------
   double CaloPulseShape::evaluate(double tDifference) const
   {
       double t = tDifference+deltaT_;
       int ibin = nSteps_ + int(t*nSteps_/digiStep_/nSteps_);

       if (ibin < 0 || ibin >= int(pulseVec_.size()-1)) return 0.0;
       double t0bin = (ibin-nSteps_)*digiStep_; //t0 is located at nSteps_
       return (pulseVec_[ibin+1]-pulseVec_[ibin])/digiStep_*(t-t0bin)+pulseVec_[ibin];
   }

   //----------------------------------------------------------------------------
   double CaloPulseShape::fromPeakToT0(double timePeak) const
   {
       return timePeak-deltaT_-0.5*digiStep_*nSteps_;
   }

   //--------------------------------
   void CaloPulseShape::diag(bool fullDiag) const
   {
       std::cout<<"Number of digi bins "<<nBinShape_<<std::endl;

       std::cout<<"Cache content "<<std::endl;
       for (auto& h : digitizedPulse_) std::cout<<h<<" ";
       std::cout<<std::endl;

       if (!fullDiag) return;

       std::cout<<"Integral val "<<std::endl;
       for (auto& h : pulseVec_) std::cout<<h<<" ";
       std::cout<<std::endl;
   }


}
