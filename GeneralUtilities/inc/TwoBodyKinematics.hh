#ifndef TwoBodyKinematics_HH
#define TwoBodyKinematics_HH
//
// Kinematics of 2 body decay.
//
// $Id: TwoBodyKinematics.hh,v 1.1 2011/05/02 16:08:34 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/02 16:08:34 $
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

#endif
