"""
Example of a python module. Contains a variable point_data,
a function add_points, and a class Point.
"""

point_data = [(1,2), (2,3), (3,4), (4,5)]

class Point:
    """
    Simple class for representing a point.
    """
    
    def __init__(self, x, y):
        """
        Create a new Point at x, y.
        """
        self.x = x
        self.y = y
        
    def translate(self, dx, dy):
        """
        Translate the point by dx and dy in the x and y directions.
        """
        self.x += dx
        self.y += dy
        
    def __add__(self, other_point):
        '''
        Adds two Point instances.
        '''
        return Point(self.x + other_point.x, 
                     self.y + other_point.y)
    
    def __str__(self):
        '''
        Return a string representation of this point instance.
        '''
        return("Point at [%f, %f]" % (self.x, self.y))

def add_points(data=None):
    """
    Example that uses point_data and the Point class.
    """
    if not data:
        data = point_data
    point_sum = None
    for point_tuple in data:
        p = Point(*point_tuple)
        if not point_sum:
            point_sum = p
        else:
            point_sum += p
    
    return point_sum