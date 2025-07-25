# Python Knot package: pythonknot

Python Knot is a Python package for creating and manipulating knots and links. It is based on the [KnotTheory](https://knot.theory.org) package for Mathematica.

## 1.Installation
only support python 3.8 or later
``` pip install pythonknot```
## 2.Features

1. Create and manipulate knots and links.
2. Compute invariants of knots and links.
3. Visualize knots and links.

## 3.Usage

    import pythonknot.alexander_poly as alexander_poly
    print(dir(alexander_poly))

    positions = alexander_poly.read_xyz("test/traj_knot31_L300_close.txt")
    print(positions.shape)
    # calculate the knot type for ring polymer
    knottype = alexander_poly.calculate_knot_type(positions, "ring")
    print(len(knottype))    

    # calculate the knot type for open polymer
    knottype = alexander_poly.calculate_knot_type(positions, "open")


    ## calculate the knot size for ring polymer
    knottype, knotsize = alexander_poly.calculate_knot_size(positions[1:100], "ring")

    print(knottype)
    print(knotsize)

    # open chain 
    positions = alexander_poly.read_xyz("test/traj_knot31_L300_open.txt")
    print(positions.shape)
    knottype, knotsize = alexander_poly.calculate_knot_size(positions, "open")
    print(knottype)
    print(knotsize)

### 3.1 HOMFLY polynomial

```python
# read traj
import pythonknot.alexander_poly as alexander_poly
positions = alexander_poly.read_xyz("traj_knot31_L300_close.txt")
print(positions.shape)

import pythonknot.homfly as homfly
crossing_number, homflr_poly = homfly.homfly_str(positions[0])  # position should be N_atom*3 array
```

### 3.2 pdb parser
load pdb trajectory, return the position of the atoms in the pdb file
```python
import pythonknot.pdb_parser as pdb_parser

import pythonknot.pdb_parser as pdb_parser

mypdb = pdb_parser.PDBParser()
mypdb.load('output.pdb')

traj = mypdb.get_all_frame_coordinates()    # with shape (N_frame, N_atom, 3)
print(traj.shape)
```

### Development

KMT algorithm, simplify the configuration of the polymer chain, return the simplified chain and the knot type of the chain

```python
#position should be a 2D array, the first dimension is the number of beads, the second dimension is the x, y, z coordinates
simplified_knot = alexander_poly.KMT_chain(positions, "open")

# if you have ring polymer, you can use the following code
simplified_knot = alexander_poly.KMT_chain(positions, "ring")
```

Extented gauss notation,
```python
# default for open chain
gauss_notation = alexander_poly.gauss_notation(simplified_knot)

# if you have ring polymer, you can use the following code, add the first point to the end of the chain
simplified_knot = np.vstack((simplified_knot, simplified_knot[0]))
gauss_notation = alexander_poly.gauss_notation(simplified_knot)

```

## TODO list
add other method of open chain knot type calculation, KMT after adding the tail
add write_xyz function to save the xyz file, modify with append mode