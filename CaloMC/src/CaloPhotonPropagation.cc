#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CaloMC/inc/CaloPhotonPropagation.hh"
#include "ConditionsService/inc/CalorimeterCalibrations.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "SeedService/inc/SeedService.hh"

#include "CLHEP/Random/RandFlat.h"
#include "TFile.h"
#include "TH2F.h"

#include <string>


namespace mu2e {

   CaloPhotonPropagation::CaloPhotonPropagation(CLHEP::HepRandomEngine& engine) : 
      timeProp_ (), 
      cdf_      (), 
      nTimeDiv_ (0), 
      dzTime_   (0), 
      randFlat_ (engine),
      lightSpeed_(300)
   {}

   //----------------------------------------------------------------------------------------------------------------------
   void CaloPhotonPropagation::buildTable()
   {
       //prepare structure for scintillating photon time propagation generation
       ConditionsHandle<CalorimeterCalibrations> calorimeterCalibrations("ignored");
       std::string fileName = calorimeterCalibrations->propagFileName();
       std::string histName = calorimeterCalibrations->propagHistName();

       TH2F *hist(0);
       TFile file(fileName.c_str());
          if (file.IsOpen()) hist = (TH2F*) file.Get(histName.c_str());
          if (!hist) throw cet::exception("CATEGORY")<<"CaloShowerStepROFromShowerStep:: Hitsogram "<<histName.c_str()
                                                     <<" from file "<<fileName.c_str()<<" does not exist\n";
	  hist->SetDirectory(0);
       file.Close();

       dzTime_   = hist->GetXaxis()->GetBinWidth(1);
       nTimeDiv_ = hist->GetNbinsY();
       for (unsigned iy=1;iy<=nTimeDiv_;++iy) timeProp_.push_back(hist->GetYaxis()->GetBinCenter(iy));
       cdf_.reserve(hist->GetNbinsX()*hist->GetNbinsY());

       for (int ix=1;ix<=hist->GetNbinsX();++ix)
       {
           float sum(1e-6);
           std::vector<float> temp;
           for (int iy=1;iy<=hist->GetNbinsY();++iy)
           {
              sum += hist->GetBinContent(ix,iy);
              temp.push_back(sum);
           }
           for (auto& val: temp) val /= sum;
           std::copy(temp.begin(),temp.end(),std::back_inserter(cdf_));
       }
       
       const Calorimeter& cal = *(GeomHandle<Calorimeter>());
       lightSpeed_            = 300.0 / cal.caloInfo().getDouble("refractiveIndex");  //in mm/ns
   }

   //----------------------------------------------------------------------------
   float CaloPhotonPropagation::propTimeSimu(float z)
   {       
       int      iz   = int(z/dzTime_);       
       float    test = randFlat_.fire(0.0,1.0);
       unsigned ibin = nTimeDiv_*iz;
       unsigned iend = ibin + nTimeDiv_;
       
       while (cdf_[ibin]<test && ibin < iend) ++ibin;
       return timeProp_[ibin-iz*nTimeDiv_];
   }

   //----------------------------------------------------------------------------
   float CaloPhotonPropagation::propTimeLine(float z)
   {       
       return z/lightSpeed_;
   }
   

}

