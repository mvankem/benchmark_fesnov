#!/usr/bin/env python3
"""Stream a .tar.gz once, reservoir-sample N files with pLDDT >= cutoff, write to a new .tar.gz."""

import io
import os
import random
import sys
import tarfile
from pathlib import Path

from tqdm import tqdm

N = 2000
SEED = 42
PLDDT_MIN = 70

input_path  = sys.argv[1]
plddt_path  = sys.argv[2]
output_path = sys.argv[3]

plddt = {}
with open(plddt_path) as f:
    for line in f:
        name, _, value = line.rstrip("\n").partition("\t")
        if value:
            plddt[name] = float(value)

rng = random.Random(SEED)
reservoir = []

with open(input_path, "rb") as raw:
    with tqdm.wrapattr(raw, "read", total=os.path.getsize(input_path),
                       unit="B", unit_scale=True, desc="reading") as fobj:
        with tarfile.open(fileobj=fobj, mode="r|gz") as tf:
            i = 0
            for member in tf:
                if not member.isfile():
                    continue
                key = Path(member.name).stem
                score = plddt.get(key)
                if score is None:
                    tqdm.write(f"warning: no pLDDT for {member.name}", file=sys.stderr)
                    continue
                if score < PLDDT_MIN:
                    continue
                if i < N:
                    reservoir.append((member, tf.extractfile(member).read()))
                else:
                    j = rng.randint(0, i)
                    if j < N:
                        reservoir[j] = (member, tf.extractfile(member).read())
                i += 1

with tarfile.open(output_path, "w:gz") as out:
    for info, data in reservoir:
        out.addfile(info, io.BytesIO(data))
