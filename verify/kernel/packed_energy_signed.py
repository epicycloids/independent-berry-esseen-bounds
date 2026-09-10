from functools import lru_cache
import json
from pathlib import Path
import time
import numpy as np
from interval_signed import SignedRectangles
from packed_energy_reference import packed_energy_reference_bounds
from complete_energy_signed import CompleteEnergyRectangles, load_complete_energy
from disk_reference_cover_engine import DiskReferenceWeights, selected, CUTOFFS
import stable_bounds

class PackedEnergyRectangles(SignedRectangles):

    def _reference(self, Llo, Lhi, dhi, taulo, tauhi):
        key = (Llo, Lhi, dhi, taulo, tauhi)
        if key in self.reference_cache:
            self.reference_cache.move_to_end(key)
            return self.reference_cache[key]
        w = self.weights
        deleted = selected(w.tlo, w.thi, 1.0, 1.0, Lhi, dhi, taulo, tauhi, deleted=True)
        E, lower, upper, eligible = ([], [], [], [])
        for j, nodes in enumerate(self.nodes):
            result = packed_energy_reference_bounds(nodes, deleted[:, j], energy_table=load_complete_energy(), L=(Llo, Lhi), d=(0.0, dhi), tau=(taulo, tauhi))
            E.append(result['E_upper'])
            lower.append(result['lower_offset'])
            upper.append(result['upper_offset'])
            eligible.append(result['eligible'])
        result = tuple((np.asarray(values).T for values in (E, lower, upper, eligible)))
        self.reference_cache[key] = result
        if len(self.reference_cache) > 16:
            self.reference_cache.popitem(last=False)
        return result
