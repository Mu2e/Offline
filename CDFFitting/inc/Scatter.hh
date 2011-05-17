#ifndef CDFFitting_Scatter_hh
#define CDFFitting_Scatter_hh
#include "CDFFitting/inc/FitAction.hh"
#include "CDFFitting/inc/Scatter.hh"
#include "CDFFitting/inc/Scatter.hh"

template <class Measureable>
class Scatter : public FitAction<Measureable> {
  
public:

  // Constructor
  Scatter();
  
  // Destructor
  virtual ~Scatter();
  
  // Access to scattering matrix
  virtual HepSymMatrix getScatterMatrix(const Measureable &  measurable) const = 0;

  virtual unsigned int getDimensionality()const =0;

  // Half of a double dispatch pair:
  virtual void applyYourselfTo(Fitter<Measureable> *theFitter) const;
  
};
#include "CDFFitting/inc/Scatter.icc"

#endif /* CDFFitting_Scatter_hh */


