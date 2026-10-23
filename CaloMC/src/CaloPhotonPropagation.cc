#include "Offline/CaloMC/inc/CaloPhotonPropagation.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "cetlib_except/exception.h"

#include "TFile.h"
#include "TH2F.h"

#include <algorithm>
#include <iterator>
#include <memory>


namespace mu2e {

   CaloPhotonPropagation::CaloPhotonPropagation(const std::string& fileName, const std::string& histName,
                                                CLHEP::HepRandomEngine& engine) :
      randFlat_(engine), fileName_(fileName), histName_(histName)
   {}


   //----------------------------------------------------------------------------------------------------------------------
   void CaloPhotonPropagation::buildTable()
   {
       constexpr float kCdfFloor = 1e-6f;
       ConfigFileLookupPolicy resolveFullPath;
       const std::string fullFileName = resolveFullPath(fileName_);

       TFile file(fullFileName.c_str());
       if (!file.IsOpen())
          throw cet::exception("CaloPhotonPropagation")
                << "cannot open propagation file " << fullFileName << "\n";

       // take ownership of the histogram so it survives the file and is freed on return
       std::unique_ptr<TH2F> hist(dynamic_cast<TH2F*>(file.Get(histName_.c_str())));
       if (!hist)
          throw cet::exception("CaloPhotonPropagation") << "histogram " << histName_
                << " not found in " << fullFileName << "\n";
       hist->SetDirectory(nullptr);
       file.Close();

       dz_       = hist->GetXaxis()->GetBinWidth(1);
       nZDiv_    = hist->GetNbinsX();
       nTimeDiv_ = hist->GetNbinsY();

       timeProp_.reserve(nTimeDiv_);
       for (unsigned iy=1; iy<=nTimeDiv_; ++iy) timeProp_.push_back(hist->GetYaxis()->GetBinCenter(iy));

       // build one normalized cumulative distribution per depth slice
       cdf_.reserve(nZDiv_*nTimeDiv_);
       for (unsigned ix=1; ix<=nZDiv_; ++ix)
       {
           std::vector<float> column;
           column.reserve(nTimeDiv_);

           float sum = kCdfFloor;
           for (unsigned iy=1; iy<=nTimeDiv_; ++iy) { sum += hist->GetBinContent(ix,iy); column.push_back(sum); }
           for (float& v : column) v /= sum;

           cdf_.insert(cdf_.end(), column.begin(), column.end());
       }

       const Calorimeter& cal = *(GeomHandle<Calorimeter>());
       lightSpeed_            = 300.0f / cal.G4Info().get<double>("refractiveIndex"); // mm/ns
   }


   //----------------------------------------------------------------------------
   float CaloPhotonPropagation::propTimeSimu(float z)
   {
       // depth slice, clamped to [0, nZDiv_-1] (negative z would be UB when cast to unsigned)
       const int      izSigned = static_cast<int>(z/dz_);
       const unsigned iz       = (izSigned <= 0) ? 0u : std::min<unsigned>(izSigned, nZDiv_-1);

       // inverse-CDF sampling: first time bin whose cumulative probability reaches the random draw
       const float test  = randFlat_.fire(0.0f,1.0f);
       const auto  first = cdf_.begin() + iz*nTimeDiv_;
       const auto  last  = first + nTimeDiv_;
       const auto  hit   = std::lower_bound(first, last, test);

       const std::size_t timeBin = std::min<std::size_t>(std::distance(first,hit), nTimeDiv_-1);
       return timeProp_[timeBin];
   }


   //----------------------------------------------------------------------------
   float CaloPhotonPropagation::propTimeLine(float z) const
   {
       return z/lightSpeed_;
   }

}
