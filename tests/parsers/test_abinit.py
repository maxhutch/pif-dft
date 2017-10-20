import unittest
from dfttopif.parsers.abinit import AbinitParser
from ..test_pif import unpack_example, delete_example
from pypif.obj.common.value import Value
import os
import shutil

class TestAbinitParser(unittest.TestCase):
        
    def get_parser(self,name):
        '''Get a PwscfParser for a certain test'''
        unpack_example(os.path.join('examples', 'abinit', name+'.tar.gz'))
        return AbinitParser(name)

    def test_Si_static(self):
        """Test that Si static example parses with AbinitParser"""
        # Parse the results
        parser = self.get_parser('abinit_Si_static')

        # Test the settings
        self.assertEquals('ABINIT', parser.get_name())

        # Test 
        self.assertEquals(0.0, parser.get_cutoff_energy().scalars[0].value)

        
if __name__ == '__main__':
    unittest.main()
