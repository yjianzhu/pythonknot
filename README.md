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
    # calculate the knot type
    knottype = alexander_poly.calculate_knot_type(positions)
    print(len(knottype))    

    knottype, knotsize = alexander_poly.calculate_knot_size(positions[1:100], "ring")

    print(knottype)
    print(knotsize)

    # open chain 
    positions = alexander_poly.read_xyz("test/traj_knot31_L300_open.txt")
    print(positions.shape)
    knottype, knotsize = alexander_poly.calculate_knot_size(positions, "open")
    print(knottype)
    print(knotsize)