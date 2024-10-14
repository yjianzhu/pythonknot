# Python Knot package: pythonknot

Python Knot is a Python package for creating and manipulating knots and links. It is based on the [KnotTheory](https://knot.theory.org) package for Mathematica.

## Installation
only support python 3.8 or later
``` pip install pythonknot```
## Features

1. Create and manipulate knots and links.
2. Compute invariants of knots and links.
3. Visualize knots and links.

## Usage

    import pythonknot.alexander_poly as alexander_poly
    print(dir(alexander_poly))

    positions = alexander_poly.read_xyz("test/traj_knot31_L300_close.txt")
    print(positions.shape)
    # calculate the knot type for ring polymer
    knottype = alexander_poly.calculate_knot_type(positions)
    print(len(knottype))    

    # calculate the knot type for open polymer
    knottype = alexander_poly.calculate_knot_type_open_chain(positions)


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
add write_xyz function to save the xyz file