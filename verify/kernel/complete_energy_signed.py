from functools import lru_cache
import json
from pathlib import Path
import time
import numpy as np
from interval_signed import SignedRectangles
from complete_energy_reference import complete_energy_reference_bounds
from energy_signed import EnergyRectangles
from disk_reference_cover_engine import DiskReferenceWeights, selected, CUTOFFS
import stable_bounds

@lru_cache(maxsize=1)
def load_complete_energy():
    path = Path(__file__).with_name('energy_table_complete.json')
    if not path.is_file():
        path = Path(__file__).parent.parent / 'farfield_phase/new_target/energy/small_frequency/energy_table_complete.json'
    result = json.loads(path.read_text())
    assert result['new_upper_interval_complete'] and result['verified_scalar_records'] == 64
    small = result['small_frequency_prefix']
    assert small['verified'] and small['right'] == '1/2' and (small['verified_rectangles'] == 25175)
    return result

class CompleteEnergyRectangles(SignedRectangles):

    def _reference(self, Llo, Lhi, dhi, taulo, tauhi):
        key = (Llo, Lhi, dhi, taulo, tauhi)
        if key in self.reference_cache:
            self.reference_cache.move_to_end(key)
            return self.reference_cache[key]
        w = self.weights
        deleted = selected(w.tlo, w.thi, 1.0, 1.0, Lhi, dhi, taulo, tauhi, deleted=True)
        E, lower, upper, eligible = ([], [], [], [])
        for j, nodes in enumerate(self.nodes):
            result = complete_energy_reference_bounds(nodes, deleted[:, j], energy_table=load_complete_energy(), L=(Llo, Lhi), d=(0.0, dhi), tau=(taulo, tauhi))
            E.append(result['E_upper'])
            lower.append(result['lower_offset'])
            upper.append(result['upper_offset'])
            eligible.append(result['eligible'])
        result = tuple((np.asarray(values).T for values in (E, lower, upper, eligible)))
        self.reference_cache[key] = result
        if len(self.reference_cache) > 16:
            self.reference_cache.popitem(last=False)
        return result
