"""Optional bounded remote replay; never invoked by the default checker."""
import json
from pathlib import Path
import modal

HERE = Path(__file__).resolve().parent if modal.is_local() else Path('/edition/verify')
ROOT = HERE.parent
app = modal.App('independent-be-certificate-replay')
image = (modal.Image.debian_slim(python_version='3.12')
    .pip_install('numpy==2.4.6', 'scipy==1.15.3', 'python-flint==0.8.0')
    .add_local_file(Path(__file__).resolve(), '/root/modal_replay.py')
    .add_local_dir(HERE/'kernel', '/edition/verify/kernel')
    .add_local_file(HERE/'scalars.py', '/edition/verify/scalars.py')
    .add_local_file(HERE/'check_saved.py', '/edition/verify/check_saved.py')
    .add_local_file(HERE/'replay.py', '/edition/verify/replay.py')
    .add_local_file(ROOT/'result.json', '/edition/result.json')
    .add_local_file(ROOT/'certificates/cover.json', '/edition/certificates/cover.json'))


@app.function(image=image, cpu=(1,1), memory=(512,1024), timeout=600, startup_timeout=60,
              max_containers=1, retries=0, scaledown_window=2, include_source=False)
def replay_short(mode: str):
    import sys
    sys.path.insert(0, '/edition/verify')
    from replay import run
    assert mode in ('topology','selected')
    return run(mode)


@app.function(image=image, cpu=(1,1), memory=(512,1024), timeout=3600, startup_timeout=60,
              max_containers=1, retries=0, scaledown_window=2, include_source=False)
def replay_one_band(index: int):
    import sys
    sys.path.insert(0, '/edition/verify')
    from replay import run
    return run('band', index)


@app.local_entrypoint()
def main(mode: str='topology', band_index: int=-1, output_dir: str='dist/replay'):
    if mode not in ('topology','selected','band'):
        raise ValueError('Choose topology, selected or band')
    if mode == 'band' and not 0 <= band_index < 767:
        raise ValueError('A band index between 0 and 766 is required')
    result = replay_one_band.remote(band_index) if mode == 'band' else replay_short.remote(mode)
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    name = f'band-{band_index:03d}.json' if mode == 'band' else f'{mode}.json'
    (output/name).write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('bands','kernel_sha256','wrapper_sha256')}, indent=2))
