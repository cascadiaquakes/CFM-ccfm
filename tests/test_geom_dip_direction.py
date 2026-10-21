"""Dip-direction validation must remain active with python -O."""
import unittest
from copy import deepcopy

from ccfm.geom import make_3d_fault_mesh


class DipDirectionTests(unittest.TestCase):
    def setUp(self):
        # An eastward trace projects south according to the right-hand rule.
        self.fault = {
            'geometry': {
                'type': 'LineString',
                'coordinates': [
                    [-120.1, 47.0], [-120.08, 47.0],
                    [-120.06, 47.0], [-120.04, 47.0],
                ],
            },
            'properties': {'dip': 45.0, 'dip_dir': 'N'},
        }

    def test_conflicting_direction_raises_value_error(self):
        with self.assertRaisesRegex(ValueError, r'projection direction .*dip direction \(N\)'):
            make_3d_fault_mesh(self.fault, check_dip_dir=True)

    def test_matching_direction_preserves_mesh(self):
        self.fault['properties']['dip_dir'] = 'S'
        unchecked = make_3d_fault_mesh(deepcopy(self.fault), check_dip_dir=False)
        checked = make_3d_fault_mesh(deepcopy(self.fault), check_dip_dir=True)
        self.assertEqual(checked, unchecked)
        self.assertGreater(len(checked), 1)
        self.assertGreater(len(checked[0]), 1)

    def test_conflicting_direction_allowed_when_check_disabled(self):
        mesh = make_3d_fault_mesh(self.fault, check_dip_dir=False)
        self.assertGreater(len(mesh), 1)
        self.assertGreater(len(mesh[0]), 1)

    def test_invalid_direction_still_rejected(self):
        self.fault['properties']['dip_dir'] = 'sideways'
        with self.assertRaisesRegex(ValueError, 'Invalid direction string'):
            make_3d_fault_mesh(self.fault, check_dip_dir=True)


if __name__ == '__main__':
    unittest.main()
