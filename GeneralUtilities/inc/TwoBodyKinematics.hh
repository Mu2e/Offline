#ifndef GeneralUtilities_TwoBodyKinematics_hh
#define GeneralUtilities_TwoBodyKinematics_hh
//
// Kinematics of 2 body decay.
//
// $Id: TwoBodyKinematics.hh,v 1.2 2011/05/17 15:41:35 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:35 $
//
// Original author Rob Kutschke

class TwoBodyKinematics{
  
public:

  TwoBodyKinematics( double m0, double m1, double m2);

  // Accept compiler generated d'tor, copy c'tor and assignment operator.

  // Accessors.
  double m0() const { return m0_;}
  double m1() const { return m1_;}
  double m2() const { return m2_;}

  double e1() const { return e1_;}
  double e2() const { return e2_;}

  double p()  const { return p_;}
  double p1() const { return p_;}
  double p2() const { return p_;}


private:

  // Copies of inputs:
  double m0_, m1_, m2_;

  // Derived quantities:
  double p_, e1_, e2_;

};

#endif /* GeneralUtilities_TwoBodyKinematics_hh */
