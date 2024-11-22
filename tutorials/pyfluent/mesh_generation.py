# from ansys.geometry.core import Point, Line, Geometry

import ansys.geometry.core
print(dir(ansys.geometry.core))


# # Define the four corner points of the rectangle
# point1 = Point(0.0, 0.0)
# point2 = Point(1.0, 0.0)
# point3 = Point(1.0, 0.5)
# point4 = Point(0.0, 0.5)

# # Create edges between consecutive points
# edge1 = Line(point1, point2)
# edge2 = Line(point2, point3)
# edge3 = Line(point3, point4)
# edge4 = Line(point4, point1)

# # Initialize the geometry object
# geom = Geometry()

# # Add the edges to the geometry
# geom.add(edge1)
# geom.add(edge2)
# geom.add(edge3)
# geom.add(edge4)

# # Optionally, save the geometry to a file
# geom.save_as('rectangle_geometry.sat')
