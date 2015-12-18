#ifndef RecoDataProducts_CaloDigi_hh
#define RecoDataProducts_CaloDigi_hh

// Original author G. Pezzullo

// C++ includes
#include <iostream>
#include <vector>

// Mu2e includes

namespace mu2e {

  struct CaloDigi{

  public:

    CaloDigi():
      _waveform(0),
      _roId    (-1),
      _nSamples(-1),
      _t0      (-1.)
    {}

    CaloDigi(const CaloDigi &CaloD):
      _waveform(CaloD.waveform()),
      _roId    (CaloD.roId()),   
      _nSamples(CaloD.nSamples()),
      _t0      (CaloD.t0())
    {}

    CaloDigi(std::vector<int> Waveform,
	     int              ROId ,
	     int              NSamples,
	     double           T0):
      _waveform(Waveform),
      _roId    (ROId),
      _nSamples(NSamples),
      _t0      (T0)
    {}

    // Accessors
    std::vector<int>      waveform()  const { return _waveform; }
    int                   roId    ()  const { return _roId;}    
    int                   nSamples()  const { return _nSamples;}
    double                t0      ()  const { return _t0;}


    // Accept compiler generated versions of d'tor, copy c'tor, assignment operator.

    // Print contents of the object.
    void print( std::ostream& ost = std::cout, bool doEndl = true ) const;

  private:

    std::vector<int>     _waveform;        //array of the samples associated with the digitezed signal
    int                  _roId;            //readout index of the sensor
    int                  _nSamples;        //total number of samples digitzed, including RO index, t0 and hit lenght
    double               _t0;              //time of the first digitezd bin of the signal

  };

 

} // namespace mu2e

#endif /* RecoDataProducts_CaloDigi_hh */
