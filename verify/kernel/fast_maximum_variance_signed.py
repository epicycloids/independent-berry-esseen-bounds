from collections import OrderedDict
import numpy as np
from flint import ctx
from fast_packed_signed import FastPackedEnergyRectangles, PreparedPrefix, PreparedGeometry, PreparedBox, fast_packed_energy_reference_bounds, frozen_reference, selected, load_complete_energy
from maximum_variance_signed import MaximumVarianceRectangles, refine_reference

class FastMaximumVarianceRectangles(FastPackedEnergyRectangles):
    score = MaximumVarianceRectangles.score

    def _reference(self, Llo, Lhi, dlo, dhi, taulo, tauhi):
        effective_precision = max(ctx.prec, 128)
        key = (Llo, Lhi, dlo, dhi, taulo, tauhi, ctx.prec, effective_precision)
        if key in self.reference_cache:
            self.reference_cache.move_to_end(key)
            return self.reference_cache[key]
        w = self.weights
        deleted = selected(w.tlo, w.thi, 1.0, 1.0, Lhi, dhi, taulo, tauhi, deleted=True)
        table = load_complete_energy()
        base = frozen_reference.reviewed_reference()
        previous_precision = ctx.prec
        E, lower, upper, eligible = ([], [], [], [])
        try:
            ctx.prec = effective_precision
            if effective_precision not in self.prepared_reference:
                self.prepared_reference[effective_precision] = (PreparedPrefix(table, base), tuple((PreparedGeometry(nodes, base) for nodes in self.nodes)))
                if len(self.prepared_reference) > 2:
                    self.prepared_reference.popitem(last=False)
            self.prepared_reference.move_to_end(effective_precision)
            prefix, geometries = self.prepared_reference[effective_precision]
            box = PreparedBox((Llo, Lhi), (dlo, dhi), (taulo, tauhi), base)
            for j, (nodes, geometry) in enumerate(zip(self.nodes, geometries)):
                result = fast_packed_energy_reference_bounds(nodes, deleted[:, j], energy_table=table, L=(Llo, Lhi), d=(dlo, dhi), tau=(taulo, tauhi), prepared_prefix=prefix, prepared_geometry=geometry, prepared_box=box)
                result = refine_reference(result, nodes, d=(dlo, dhi), tau=(taulo, tauhi))
                E.append(result['E_upper'])
                lower.append(result['lower_offset'])
                upper.append(result['upper_offset'])
                eligible.append(result['eligible'])
        finally:
            ctx.prec = previous_precision
        result = tuple((np.asarray(values).T for values in (E, lower, upper, eligible)))
        self.reference_cache[key] = result
        if len(self.reference_cache) > 16:
            self.reference_cache.popitem(last=False)
        return result
