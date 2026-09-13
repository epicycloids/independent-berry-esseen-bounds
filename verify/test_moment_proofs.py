"""Controls for complete nested moment covers and fixed signed consumers."""
from copy import deepcopy
from fractions import Fraction as Q
import unittest
import zipfile
from check_saved import ROOT,load_inputs,load_band
from moment_proofs import validate,signed_edges


class MomentProofTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cover,_=load_inputs();cls.samples={}
        def visit(proof,lo,hi):
            rule=proof['rule'];cls.samples.setdefault(rule,(proof,lo,hi))
            if rule=='evaluation':
                recipe=proof['recipe'];cls.samples.setdefault(recipe['method'],(proof,lo,hi))
                if recipe.get('consumer')=='cantelli':cls.samples.setdefault('cantelli',(proof,lo,hi))
            if rule=='cover':
                for child in proof['regions']:visit(child['proof'],lo,hi)
            elif rule=='feasible_b':visit(proof['proof'],lo,hi)
        wanted={'cover','empty','signed','disk_signed','cantelli'}
        with zipfile.ZipFile(ROOT/cover['band_archive']) as archive:
            for row in cover['bands'][:144]:
                band=load_band(archive,row)
                for leaf in band['leaves']:
                    if 'proof' in leaf:visit(leaf['proof'],band['Llo'],band['Lhi'])
                if wanted<=cls.samples.keys():break
        if not wanted<=cls.samples.keys():raise AssertionError('Missing actual proof cases')

    def changed(self,key,edit):
        source,lo,hi=self.samples[key];proof=deepcopy(source);edit(proof)
        with self.assertRaises((ValueError,AssertionError,KeyError)):validate(proof,lo,hi)

    def test_actual_complete_samples(self):
        for proof,lo,hi in self.samples.values():validate(proof,lo,hi)

    def test_missing_closed_child(self):
        self.changed('cover',lambda p:p['regions'].pop())

    def test_wrong_partition_maximum(self):
        self.changed('cover',lambda p:p.update(upper=str(Q(p['upper'])+1)))

    def test_changed_inherited_coordinate(self):
        def edit(p):p['regions'][0]['box'][4]+=0.001
        self.changed('cover',edit)

    def test_missing_threshold_interval(self):
        self.changed('signed',lambda p:p['recipe']['theta_edges'].pop(5))

    def test_missing_outside_tail_consumer(self):
        self.changed('cantelli',lambda p:p['recipe'].pop('consumer'))

    def test_missing_disk_frequency_cell(self):
        self.changed('disk_signed',lambda p:p['recipe']['contexts'][0]['coarse_cells'].pop())

    def test_changed_disk_signed_normalization(self):
        def edit(p):p['recipe']['contexts'][0]['coarse_cells'][0]['joint_upper']='0'
        self.changed('disk_signed',edit)

    def test_missing_refined_rectangle(self):
        def edit(p):
            c=next(c for c in p['recipe']['contexts'] if c['threshold_pieces'])
            c['threshold_pieces'][0]['cells'][0]['pieces'].pop()
        self.changed('disk_signed',edit)

    def test_resource_exclusion_needs_strict_gap(self):
        proof=dict(rule='empty',box=[0.,1.,0.,1.,0.,1.],upper='0',condition='remainder_resource')
        with self.assertRaises(ValueError):validate(proof,.5,1.)

    def test_third_moment_exclusion_needs_strict_gap(self):
        proof=dict(rule='empty',box=[.5,1.,0.,1.,0.,1.],upper='0',condition='b_squared')
        with self.assertRaises(ValueError):validate(proof,.5,1.)

    def test_threshold_covers_keep_full_outside_endpoints(self):
        for count,end in ((129,Q(2)),(139,Q(69,32)),(145,Q(69,32))):
            edges=signed_edges(count);self.assertEqual((len(edges),edges[0],edges[-1]),(count,-end,end))


if __name__=='__main__':unittest.main()
