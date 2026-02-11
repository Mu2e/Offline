#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CaloMC/inc/CaloPhotonPropagation.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/ConditionsService/inc/ConditionsHandle.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/SeedService/inc/SeedService.hh"

#include "CLHEP/Random/RandFlat.h"
#include "TFile.h"
#include "TH2F.h"

#include <string>


namespace mu2e {

   CaloPhotonPropagation::CaloPhotonPropagation(const std::string& fileName, const std::string& histName, CLHEP::HepRandomEngine& engine) :
      timeProp_ (),
      cdf_      (),
      nTimeDiv_ (0),
      nZDiv_ (0),
      dzTime_   (0),
      randFlat_ (engine),
      fileName_(fileName),
      histName_(histName),
      lightSpeed_(300)
   {}

   //----------------------------------------------------------------------------------------------------------------------
   void CaloPhotonPropagation::buildTable()
   {
       ConfigFileLookupPolicy resolveFullPath;
       std::string fullFileName = resolveFullPath(fileName_);

       TH2F *hist(0);
       TFile file(fullFileName.c_str());
          if (file.IsOpen()) hist = (TH2F*) file.Get(histName_.c_str());
          if (!hist) throw cet::exception("CATEGORY")<<"CaloROStepMaker:: Hitsogram "<<histName_.c_str()
                                                     <<" from file "<<fileName_.c_str()<<" does not exist\n";
          hist->SetDirectory(0);
       file.Close();

       dzTime_   = hist->GetXaxis()->GetBinWidth(1);
       nTimeDiv_ = hist->GetNbinsY();
       nZDiv_ = hist->GetNbinsX();
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
       unsigned iz   = z/dzTime_;
       if(iz>=nZDiv_) iz = nZDiv_ - 1;
       float    test = randFlat_.fire(0.0,1.0);
       unsigned ibin = nTimeDiv_*iz;
       unsigned iend = ibin + nTimeDiv_ - 1;

       while (cdf_[ibin]<test && ibin < iend) ++ibin;
       return timeProp_[ibin-iz*nTimeDiv_];
   }

   //----------------------------------------------------------------------------
   float CaloPhotonPropagation::propTimeLine(float z)
   {
       return z/lightSpeed_;
   }


}

