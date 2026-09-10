from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
R27 = HERE.parent if (HERE.parent / 'fast_maximum_cover_engine.py').is_file() else HERE
for folder in (HERE, R27, R27.parent / '26', R27 / 'dual_logarithm', R27 / 'disk_table', R27 / 'reference_ode', R27 / 'reference_ode/efficiency', R27 / 'farfield_phase/new_target', R27 / 'farfield_phase/new_target/energy', R27 / 'farfield_phase/new_target/energy/small_frequency'):
    if folder.is_dir() and str(folder) not in sys.path:
        sys.path.insert(0, str(folder))
