# Python Knot package: pythonknot

Python Knot is a package for creating and manipulating knots and links.

## 1. Installation

Python 3.8+ is required.

```bash
pip install pythonknot
```

### 1.1 Alexander backend (Rust + PyO3)

`pythonknot.alexander_poly` uses the Rust native backend.

- packaged table: installed with wheel at `pythonknot/data/table_knot_Alexander_polynomial.txt` and auto-used on import
- legacy fallback path: `~/.local/share/rust_knot/table_knot_Alexander_polynomial.txt`
- override table: `PYTHONKNOT_ALEXANDER_TABLE=/path/to/table.txt`

`homfly` C/C++ extension is removed from build; wheels now contain only the Rust backend.

## 2. API Summary (`pythonknot.alexander_poly`)

- `read_xyz(filename) -> np.ndarray[F, N, 3]`
- `read_pdb(filename, atom_filter="all") -> np.ndarray[F, N, 3]`
- `write_xyz(filename, input_data, append=False) -> None`
- `knot_type(input_data, chain_type="ring", threads=None) -> List[str]`
- `knot_size(input_data, chain_type="ring", threads=None) -> Tuple[List[str], List[List[int]]]`
- `kmt(input_data, chain_type="ring") -> np.ndarray[N2, 3]`

Notes:
- `input_data` can be XYZ path (`str`), `(N,3)` array, or `(F,N,3)` array for `knot_type` / `knot_size`.
- `read_pdb(..., atom_filter="all")` supports `"all"` and `"ca"`.
- `kmt` requires a `(N,3)` array.
- `threads=None` means single-thread path; `threads>1` enables rayon parallel processing.
- `knot_type` and `knot_size` release Python GIL during heavy Rust computation.

## 3. Usage

```python
import numpy as np
import pythonknot.alexander_poly as alexander_poly

positions = alexander_poly.read_xyz("test/traj_knot31_L300_close.txt")
print(positions.shape)

# knot type
types_ring = alexander_poly.knot_type(positions, chain_type="ring", threads=8)
types_open = alexander_poly.knot_type(positions, chain_type="open", threads=8)
print(types_ring[:5])

# knot size
types, sizes = alexander_poly.knot_size(positions[:100], chain_type="ring", threads=8)
print(types[:3], sizes[:3])

# KMT simplification for one frame (N,3)
simplified = alexander_poly.kmt(positions[0], chain_type="ring")
print(simplified.shape)

# write (overwrite by default)
alexander_poly.write_xyz("out.xyz", positions, append=False)

# append frames
alexander_poly.write_xyz("out.xyz", positions[:10], append=True)

# read PDB trajectory directly
pdb_all = alexander_poly.read_pdb("output.pdb", atom_filter="all")
pdb_ca = alexander_poly.read_pdb("output.pdb", atom_filter="ca")
print(pdb_all.shape, pdb_ca.shape)
```

## 4. PDB reading

Use `pythonknot.alexander_poly.read_pdb(...)` directly.

## 5. Cross-platform wheel build

Build is platform-native and Rust-only.

Linux (manylinux, cp38-cp314):
```bash
./compile_chain.sh
```

macOS arm64 (run on macOS arm64 host):
```bash
for py in python3.8 python3.9 python3.10 python3.11 python3.12 python3.13 python3.14; do
  "$py" -m pip install -U pip setuptools wheel
  "$py" setup.py bdist_wheel
done
```

Windows x64 (PowerShell, run on Windows host):
```powershell
$versions = @("3.8","3.9","3.10","3.11","3.12","3.13","3.14")
foreach ($v in $versions) {
  py -$v -m pip install -U pip setuptools wheel
  py -$v setup.py bdist_wheel
}
```
