#ifndef GeneralUtilities_TwoBodyKinematics_hh
#define GeneralUtilities_TwoBodyKinematics_hh
//
// Kinematics of 2 body decay.
//
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
