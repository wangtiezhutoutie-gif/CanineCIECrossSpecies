import json
import sys
from pathlib import Path

import h5py


def describe(path: str) -> None:
    p = Path(path)
    print(f"FILE\t{p}\t{p.stat().st_size}")
    with h5py.File(p, "r") as f:
        print("ROOT_ATTRS\t" + json.dumps({k: str(v) for k, v in f.attrs.items()}, sort_keys=True))
        for key in f.keys():
            obj = f[key]
            shape = getattr(obj, "shape", None)
            attrs = {k: str(v) for k, v in obj.attrs.items()}
            print(f"ROOT\t{key}\t{type(obj).__name__}\t{shape}\t{json.dumps(attrs, sort_keys=True)}")
        for group_name in ("obs", "var", "layers", "raw", "obsm", "uns"):
            if group_name not in f:
                continue
            group = f[group_name]
            print(f"GROUP\t{group_name}\t{list(group.keys())}")
            for key in group.keys():
                obj = group[key]
                shape = getattr(obj, "shape", None)
                attrs = {k: str(v) for k, v in obj.attrs.items()}
                children = list(obj.keys()) if isinstance(obj, h5py.Group) else None
                print(f"ITEM\t{group_name}/{key}\t{type(obj).__name__}\t{shape}\t{children}\t{json.dumps(attrs, sort_keys=True)}")


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        describe(arg)
