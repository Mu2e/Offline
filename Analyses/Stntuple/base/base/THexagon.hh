///////////////////////////////////////////////////////////////////////////////
// HexSize - distance between the 2 flat sides
///////////////////////////////////////////////////////////////////////////////
#ifndef test_hexmap_hh
#define test_hexmap_hh

// #include "TCanvas.h"
#include "TPolyLine.h"
#include "TMath.h"
#include "TVector2.h"
//-----------------------------------------------------------------------------
class THexagon: public TObject {
public:
  double fSize;
  double fX0;
  double fY0;
  
  int    fLineColor;
  int    fFillColor;
  int    fFillStyle;
//-----------------------------------------------------------------------------
// methods
//-----------------------------------------------------------------------------
  THexagon();
  ~THexagon();
  THexagon(double HexSize, double X = 0, double Y = 0);

  double X0    () const { return fX0; }
  double Y0    () const { return fY0; }
  double Radius() const { return sqrt (fX0*fX0+fY0*fY0); }

  double Dist(double X, double Y) const { return sqrt((fX0-X)*(fX0-X)+(fY0-Y)*(fY0-Y)); }

  void   SetSize(double HexSize) { fSize = HexSize; }
  void   SetPos (double X, double Y) { fX0 = X; fY0 = Y; }


  void  SetLineColor(int Color) { fLineColor = Color; }
  void  SetFillColor(int Color) { fFillColor = Color; }
  void  SetFillStyle(int Style) { fFillStyle = Style; }

  void  Paint(Option_t* Opt="") ;

  static int Test(double HexSize = 30., double RIn=360., double ROut=670.);

  ClassDef(THexagon,0)

};


#endif
