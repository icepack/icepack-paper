import argparse
import numpy as np
import pygmsh

parser = argparse.ArgumentParser()
parser.add_argument('--output')
args = parser.parse_args()

R = 200e3
δx = 5e3
geometry = pygmsh.built_in.Geometry()

x1 = geometry.add_point([-R, 0, 0], lcar=δx)
x2 = geometry.add_point([+R, 0, 0], lcar=δx)

center1 = geometry.add_point([0, 0, 0,], lcar=δx)
center2 = geometry.add_point([0, -4 * R, 0], lcar=δx)

arcs = [geometry.add_circle_arc(x1, center1, x2),
        geometry.add_circle_arc(x2, center2, x1)]

line_loop = geometry.add_line_loop(arcs)
plane_surface = geometry.add_plane_surface(line_loop)

physical_lines = [geometry.add_physical(arc) for arc in arcs]
physical_surface = geometry.add_physical(plane_surface)

with open(args.output, 'w') as geo_file:
    geo_file.write(geometry.get_code())
