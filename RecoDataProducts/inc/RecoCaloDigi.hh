#ifndef RecoDataProducts_RecoCaloDigi_hh
#define RecoDataProducts_RecoCaloDigi_hh

// $Id: $
// $Author:  $
// $Date:  $
//

// C++ includes
#include <iostream>
#include <vector>

#include "RecoDataProducts/inc/CaloDigiCollection.hh"

namespace mu2e {

  class RecoCaloDigi {

  public:

    RecoCaloDigi():
      _ROid    (-1), 
      _caloDigi(CaloDigi()), 
      _edep    (-1),
      _amplitude(-1),
      _time    (-1),
      _tChi2   (-1),
      _psd     (-1) {} 

    RecoCaloDigi(const RecoCaloDigi &RecoCaloD):
      _ROid     (RecoCaloD.ROid()), 
      _caloDigi (RecoCaloD.RODigi()), 
      _edep     (RecoCaloD.edep()),
      _amplitude(RecoCaloD.amplitude()),
      _time     (RecoCaloD.time()),
      _tChi2    (RecoCaloD.tChi2()),
      _psd      (RecoCaloD.psd()) {} 

    RecoCaloDigi(CaloDigi CaloD,
		 double   Edep, 
		 double   Amplitude, 
		 double   Time, 
		 double   TChi2, 
		 double   Psd):
      _ROid    (CaloD.roId()), 
      _caloDigi(CaloD), 
      _edep    (Edep),
      _amplitude(Amplitude),
      _time    (Time),
      _tChi2   (TChi2),
      _psd(Psd) {}


    //Accessors
    int       ROid     () const { return _ROid;}               //photosensor id      
    CaloDigi  RODigi   () const { return _caloDigi;}           // starting position of the hit in the DigiOutput vector 
       
    double    edep     () const { return _edep;}               // reconstructed energy in MeV
    double    amplitude() const { return _amplitude;}          //  waveform amplitude in mV
    double    time     () const { return _time;}               // reconstructed time in ns
    double    tChi2    () const { return _tChi2;}              // chi2 on time fit
       
    double    psd() const { return _psd;}     // hit from pileup event flag
    	
    
  private:
    int       _ROid;
    CaloDigi  _caloDigi;          // starting position of the hit in the DigiOutput vector 
       
    double    _edep;               // reconstructed energy
    double    _amplitude;          // waveform amplitude
    double    _time;               // reconstructed time
    double    _tChi2;               // chi2 on time fit
       
    double    _psd;     // hit from pileup event flag
  };
}

#endif /* RecoDataProducts_RecoCaloDigi_hh */
