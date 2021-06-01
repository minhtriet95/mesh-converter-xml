# Copyright (C) 2021 by Triet Pham
#
# February 10, 2021

from dolfin import *
from mshr import *
import os
import timeit
import meshio
import numpy as np

start_time = timeit.default_timer()

# -------------------- Extract file extension -------------------- #
filePath = 'cellAffine.msh'
fileName, fileExtension = os.path.splitext(filePath)

# -------------------- Read mesh -------------------- #
msh = meshio.read(filePath)

dim = None
element_count = dict()  # Store element info into dict()
element_types = list()  # Element types
for cell in msh.cells:
    element_count[str(cell.type)] = len(cell.data)
    
    #! Only support 'line', 'triangle', and 'tetra'
    element_types.append(cell.type)
    dim = 3 if cell.type == 'tetra' else 2

points = msh.points
if dim == 2: points = np.delete(points, 2, axis=1)  # Remove z column

dash = '-' * 40
print(dash)
print('{:^40}'.format('MESH INFO'))
print(dash)

print('{:<25}{:>15}'.format('Case:', str(dim) + '-D'))
print('{:<25}{:>15}'.format('Number of nodes:', len(points)))
print('{:<25}'.format('Number of elements:'))
for ele_type, ele_count in element_count.items():
    print('{:<15}{:<15}{:>10}'.format('', ele_type, ele_count))

print(dash)


# -------------------- Convert to xdmf files -------------------- #
meshio.write("xdmf/mesh.xdmf", meshio.Mesh(
    points, cells={element_types[1]: msh.cells_dict[element_types[1]]}))
meshio.write("xdmf/subdomains.xdmf", meshio.Mesh(
    points, cells={element_types[1]: msh.cells_dict[element_types[1]]},
    cell_data={"name_to_read": [msh.cell_data["gmsh:physical"][1]]}))
meshio.write("xdmf/boundaries.xdmf", meshio.Mesh(
    points, cells={element_types[0]: msh.cells_dict[element_types[0]]},
    cell_data={"name_to_read": [msh.cell_data["gmsh:physical"][0]]}))


# -------------------- Reload xdmf files -------------------- #
mesh = Mesh()
with XDMFFile("xdmf/mesh.xdmf") as infile:
    infile.read(mesh)

# Physical (subdomain) regions
mvc_block = MeshValueCollection("size_t", mesh, dim)
with XDMFFile("xdmf/subdomains.xdmf") as infile:
    infile.read(mvc_block, "name_to_read")
subdomains = MeshFunction("size_t", mesh, mvc_block)

# Facet (boundary) regions
mvc_facet = MeshValueCollection("size_t", mesh, dim - 1)
with XDMFFile("xdmf/boundaries.xdmf") as infile:
    infile.read(mvc_facet, "name_to_read")
boundaries = MeshFunction("size_t", mesh, mvc_facet)

cell_data = msh.cell_data["gmsh:physical"][1]
for index, value in enumerate(mesh.cells()):
    mesh.domains().set_marker([index, cell_data[index]], dim)

# Add sideset info to boundaries
sideset_list = set(msh.cell_data["gmsh:physical"][0])
cell_index_list = list()

for sideset_id in sideset_list:
    cell_index_list.append(boundaries.where_equal(sideset_id))

boundaries.set_all(0)
    
for index, sideset_id in enumerate(sideset_list):
    for cell_index in cell_index_list[index]:
        boundaries.set_value(cell_index, sideset_id)

# sideset_1 = boundaries.where_equal(1)
# sideset_2 = boundaries.where_equal(2)
# sideset_3 = boundaries.where_equal(3)

# boundaries.set_all(0)

# for cell_index in sideset_1:
#     boundaries.set_value(cell_index, 1)
# for cell_index in sideset_2:
#     boundaries.set_value(cell_index, 2)
# for cell_index in sideset_3:
#     boundaries.set_value(cell_index, 3)


# -------------------- Export to .xml files -------------------- #
File(fileName + ".xml") << mesh
File(fileName + "_physical_region.xml") << subdomains
File(fileName + "_facet_region.xml") << boundaries

stop_time = timeit.default_timer()
run_time = round(stop_time - start_time, 2)

print('{:^40}'.format(f'FINISH CONVERT IN {run_time}s'))
print(dash)
