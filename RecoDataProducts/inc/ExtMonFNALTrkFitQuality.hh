// This class keeps info pertinent to individual track fit, such as
// chi2.  Quantities that characterize a track in an ensemble (such as
// number of shared hits for a track) do not belong here.
//
// Original author Andrei Gaponenko
//
//
// $Id: ExtMonFNALTrkFitQuality.hh,v 1.2 2012/11/01 23:43:39 gandr Exp $
// $Author: gandr $
// $Date: 2012/11/01 23:43:39 $

#ifndef RecoDataProducts_ExtMonFNALTrkFitQuality_hh
#define RecoDataProducts_ExtMonFNALTrkFitQuality_hh

#include <ostream>

namespace mu2e {

  //================================================================
  class ExtMonFNALTrkFitQuality {
  public:

    double chi2() const { return chi2_;}
    int    ndf()  const { return ndf_;}

    ExtMonFNALTrkFitQuality(int ndf, double chi2)
      : chi2_(chi2), ndf_(ndf)
    {}

    // ctr for genreflex persistency
    ExtMonFNALTrkFitQuality()
      : chi2_(), ndf_()
    {}

  private:
    double chi2_;
    int ndf_;
  };


  std::ostream& operator<<(std::ostream&os, const ExtMonFNALTrkFitQuality& q);

} // namespace mu2e

#endif /* RecoDataProducts_ExtMonFNALTrkFitQuality_hh */
