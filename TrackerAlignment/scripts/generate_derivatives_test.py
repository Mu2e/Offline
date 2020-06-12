import unittest 

import generate_derivatives

class NumericalTests(unittest.TestCase):
    def test_DOCA(self):
        """
        Test that the implemented DOCA function returns the expected value for a given case
        """

        pass

    def test_alignment_function(self):
        """
        Test that the straw position and direction alignment behaves in line with the Mu2e Offline code (AlignedTrackerMaker)
        within some tolerance

        """

        aligned_wpos = aligned_wpos.subs({
            'dx': 1.50269,
            'dy': -1.94275,
            'dz': 3.27363,

            'a': 0,
            'b': 0,
            'g': 0,

            'panel_dx': 0,
            'panel_dy': 0,
            'panel_dz': 0,
            'panel_a': 0,
            'panel_b': 0,
            'panel_g': 0,

            'wx': 367.052,  # wire pos before (straw 0 plane 0 panel 0)
            'wy': 98.3512,
            'wz': -1490.37,

            'panel_x': 368.561,  # panel straw0mp
            'panel_y': 98.7556,
            'panel_z': -1493.08,

            'ppx': 0,  # plane origin
            'ppy': 0,
            'ppz': -1504.35
        })

        
        pass

    def test_derivatives(self):
        """
        Test that the analytical solutions for the derivatives are consistent with numerical
        results (TODO)

        """

        pass