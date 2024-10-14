import pythonknot
import numpy
print(dir(pythonknot))  
import pythonknot.xyz as xyz
print(dir(xyz))

import pythonknot.alexander_poly as alexander_poly
print(dir(alexander_poly))

#alexander_poly.calculate_knot_type("traj_knot31_L300_close.txt")

# positions = alexander_poly.read_xyz("traj_knot31_L300_close.txt")

# print(positions.shape)
# simplified_knot = alexander_poly.KMT_chain(positions[0], "ring")

# print(simplified_knot)
# simplified_knot = numpy.expand_dims(simplified_knot, axis=0)
# print(alexander_poly.calculate_knot_type(simplified_knot))
# knottype = alexander_poly.calculate_knot_type(positions)
# print(len(knottype))    

# knottype, knotsize = alexander_poly.calculate_knot_size(positions[1:100], "ring")

# print(knottype)
# print(knotsize)

positions = alexander_poly.read_xyz("traj_knot31_L300_open.txt")

print(positions.shape)
# knottype, knotsize = alexander_poly.calculate_knot_size(positions, "open")
# print(knottype)
# print(knotsize)
simplified_knot = alexander_poly.KMT_chain(positions[0], "open")
print(simplified_knot)

# 扩充维度

# simplified_knot = numpy.expand_dims(simplified_knot, axis=0)
# print(alexander_poly.calculate_knot_type_open_chain(simplified_knot))

gauss_notation = alexander_poly.gauss_notation(simplified_knot)
print(gauss_notation)
