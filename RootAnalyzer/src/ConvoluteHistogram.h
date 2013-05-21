// header for convoluting histogram code

// $Id: ConvoluteHistogram.h,v 1.1 2013/05/21 14:25:44 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2013/05/21 14:25:44 $
//
// Original author Kyle Knoepfel
//

class TH1;
class ConvoluteHistogram {

 public:

  static TH1* GetConvHistogram( const TH1& histo, const TH1& hConv );
  static TH1* GetConvHistogram( const TH1& histo, 
                                const std::string functionLoc,
                                const double halfwidth );
  // use compiler generated c'stor, d'stor

 private:

  static const int _nSamplings = 10000;

};
#ifdef __CINT__
#pragma link C++ class ConvoluteHistogram+;
#endif
