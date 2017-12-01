from lens import Lens
import unittest

class TestLens(unittest.TestCase):

	def setUp(self):
		self.lens1 = Lens(123456789,30.0,60.0,100.0,100.0,2009.0)
		self.lens2 = Lens(11111111,31.0,61.0,101.0,102.0,2009.1)

	def test_id_init(self):
		self.assertEqual(self.lens1._id,123456789)
		self.assertEqual(self.lens2._id,11111111)
	
	def test_ra_init(self):
		self.assertEqual(self.lens1._ra_0,30.0)
		self.assertEqual(self.lens2._ra_0,31.0)

	def test_dec_init(self):
		self.assertEqual(self.lens1._dec_0,60.0)
		self.assertEqual(self.lens2._dec_0,61.0)

	def test_pmra_init(self):
		self.assertEqual(self.lens1._pmra,100.0)
		self.assertEqual(self.lens2._pmra,101.0)
		
	def test_pmdec_init(self):
		self.assertEqual(self.lens1._pmdec,100.0)
		self.assertEqual(self.lens2._pmdec,102.0)
		
	def test_epoch_init(self):
		self.assertEqual(self.lens1._epoch_0,2009.0)
		self.assertEqual(self.lens2._epoch_0,2009.1)

	def test_get_eq_coord_epoch(self):
		self.assertAlmostEqual(self.lens1.get_eq_coords_at_epoch(2025.0)[0],15.000444444)
		self.assertAlmostEqual(self.lens1.get_eq_coords_at_epoch(2025.0)[1],60.000444444)
		self.assertAlmostEqual(self.lens1.get_eq_coords_at_epoch(1956.0)[0],14.99852778)
		self.assertAlmostEqual(self.lens1.get_eq_coords_at_epoch(1956.0)[1],59.99852779)
		self.assertAlmostEqual(self.lens1.get_eq_coords_at_epoch(2017.0)[0],15.00022222)
		self.assertAlmostEqual(self.lens1.get_eq_coords_at_epoch(2017.0)[1],60.00022222)

	#def test_get_seprations_at_epoch(self):
		self.assertAlmostEqual(self.lens1.get_angular_separation_at_epoch(2017.0,30.00069760,60.00024),90.5,places=1)
	
		
if __name__ == '__main__':
	suite = unittest.TestLoader().loadTestsFromTestCase(TestLens)
	unittest.TextTestRunner(verbosity=2).run(suite)

