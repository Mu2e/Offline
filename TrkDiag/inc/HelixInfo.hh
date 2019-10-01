//
//  Struct for HelixSeed info
//
#ifndef TrkDiag_HelixInfo_hh
#define TrkDiag_HelixInfo_hh
namespace mu2e {
  struct HelixInfo {
    int _nch, _ncha;
    int _nsh, _nsha;
    unsigned _flag;
    float _t0err;
    float _mom;
    float _chi2xy;
    float _chi2fz;
    float _ecalo;
    void reset() { _nch = _nsh = _ncha = _nsha = -1; _flag = 0;
      _t0err = _mom = _chi2xy = _chi2fz = _ecalo =-1.0; }
    static std::string leafnames() { 
      static std::string leaves
	= std::string("nch/I:ncha/I:nsh/I:nsha/I:flag/i:t0err/F:mom/F:chi2xy/F:chi2fz/F:ecalo/F");
      return leaves;
    }
  };
}
#endif
