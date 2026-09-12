"""Explicit, pinned source acquisition. Never import or execute upstream code."""
from __future__ import annotations
from pathlib import Path
import subprocess
from .provenance import freeze_checkout,atomic_json,git
from .registry import import_predecessor_registry

SOURCE_COMMIT='3ed32805b38ce34216b34888f6539e3538e90fb9'
SOURCE_URL='https://github.com/uskutsav-cpu/selfdual-5form-invariants.git'
REQUIRED_FILES=('results/10d_order8.json','results/10d_order10.json','results/10d_order12.json',
    'src/sdinv/stress.py','src/sdinv/interaction.py','src/sdinv/invariant_registry.py',
    'docs/interacting_stress_tensor.md','docs/assumptions_limitations_and_open_questions.md')


def acquire_source(checkout,*,fetch=False):
    checkout=Path(checkout).resolve()
    if not checkout.exists():
        if not fetch:
            raise FileNotFoundError('source checkout is missing; --fetch explicitly authorizes cloning the pinned public predecessor')
        checkout.parent.mkdir(parents=True,exist_ok=True)
        subprocess.run(['git','clone','--no-checkout',SOURCE_URL,str(checkout)],check=True)
        subprocess.run(['git','-C',str(checkout),'checkout','--detach',SOURCE_COMMIT],check=True)
    # Existing source worktrees are never reset or overwritten.
    manifest=freeze_checkout(checkout,REQUIRED_FILES,SOURCE_COMMIT)
    return checkout,manifest


def export_source(checkout,output,fetch=False):
    root,manifest=acquire_source(checkout,fetch=fetch)
    registry=import_predecessor_registry(root,12)
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    atomic_json(output/'source_manifest.json',manifest)
    atomic_json(output/'registry.json',registry.to_json())
    return {'manifest':str(output/'source_manifest.json'),'registry':str(output/'registry.json'),
            'source_commit':SOURCE_COMMIT,'status':'source_frozen_not_physics_recomputed'}
