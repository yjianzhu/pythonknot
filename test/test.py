import pythonknot
print(dir(pythonknot))  
import pythonknot.xyz as xyz
print(dir(xyz))

import pythonknot.alexander_poly as alexander_poly
print(dir(alexander_poly))

#alexander_poly.calculate_knot_type("test/traj_knot31_L300_close.txt")

positions = alexander_poly.read_xyz("test/traj_knot31_L300_close.txt")

print(positions.shape)
knottype = alexander_poly.calculate_knot_type(positions)
print(len(knottype))    

knottype, knotsize = alexander_poly.calculate_knot_size(positions[1:100], "ring")

print(knottype)
print(knotsize)

positions = alexander_poly.read_xyz("test/traj_knot31_L300_open.txt")

print(positions.shape)
knottype, knotsize = alexander_poly.calculate_knot_size(positions, "open")
print(knottype)
print(knotsize)