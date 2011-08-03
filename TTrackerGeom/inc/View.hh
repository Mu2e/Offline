#ifndef TrackerGeom_View_hh
#define TrackerGeom_View_hh

//
// Holds the discrete angular view number of a panel 
//
// $Id: View.hh,v 1.1 2011/08/03 18:31:25 mf Exp $
// $Author: mf $
// $Date: 2011/08/03 18:31:25 $
//
// Original author Mark Fischler
//
// Our definitions:  The "angle" (or phi) of a panel is determined as 
// an angle from Y of a line perpendicular to the wire directions.  
// Let us be more specific:
//
// The direction of a wire is the direction given by Straw::getDirection()
// for the straw of that wire.  We believe this is the direction TOWARDS
// (rather than away from) the readout end of the wire; the end at which t is 
// recorded.  (The other end pulse time is t+deltaT.)
//
// Then phi for a panel is 90 degrees less than the direction of its wires.
//
// The View of a panel is in principle ((phi - phi0) mod 360 degrees)) / 30; 
// this ranges  from 0 to 11.  We will call this the "hour" of the View.
// Here phi0 is some reference orientation dictated by the rotational 
// orientations of the stations.  Our definition is that the lowest-Z face
// of each station has one of its panels at a an angle of phi0 and a view of 0.
//
// There is a rounding issue here.  If we take the floor of 
// ((phi - phi0) mod 360 degrees)) / 30 as the view number, then for 
// phi0 = 0 each of the angles is an exact multiple of 30 degrees, thus
// susceptible to slight round-down causing a wrong view assignment.
// If we take the nearest integer, then for phi0 = 15 degrees a similer
// problem arises.  Currently, the stations are planned for 15 degrees;
// they might be at 0 degrees.  So neither scheme works well.  
//
// Our compromise is that view = floor(((phi - phi0 + 5) mod 360 degrees)) / 30.   

namespace mu2e {

  typedef int View;

}

#endif // TrackerGeom_View_hh
