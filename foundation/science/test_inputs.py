"""Adversarial checks of the new standalone input boundary and provenance."""
import ast,copy,hashlib,importlib.util,json,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('standalone',ROOT/'rank81/verify.py')
v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
class InputTests(unittest.TestCase):
    def setUp(self):self.data=json.loads((ROOT/'rank81/input.json').read_text())
    def test_valid(self):v.validate_input(self.data)
    def test_changed_metric(self):
        self.data['convention']['metric'][0]=1
        with self.assertRaises(AssertionError):v.validate_input(self.data)
    def test_composite(self):
        self.data['points'][0]['prime']=65535
        with self.assertRaises(AssertionError):v.validate_input(self.data)
    def test_missing_graph(self):
        self.data['graphs'].pop()
        with self.assertRaises(AssertionError):v.validate_input(self.data)
    def test_loop(self):
        self.data['graphs'][0]['graph']['edges'][0][1]=0
        with self.assertRaises(AssertionError):v.validate_input(self.data)
    def test_bad_coordinates(self):
        self.data['points'][0]['coordinates'][0]=-1
        with self.assertRaises(AssertionError):v.validate_input(self.data)
    def test_provenance(self):
        source=(ROOT/'rank81/verify.py').read_text();tree=ast.parse(source)
        bodies={n.name:ast.get_source_segment(source,n) for n in tree.body if isinstance(n,(ast.FunctionDef,ast.ClassDef))}
        for record in json.loads((ROOT/'rank81/implementation_provenance.json').read_text())['symbols']:
            self.assertEqual(hashlib.sha256(bodies[record['symbol']].encode()).hexdigest(),record['source_sha256'])
if __name__=='__main__':unittest.main()
