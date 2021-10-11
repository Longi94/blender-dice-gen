import math
import bpy
import os
from typing import List
from math import sqrt, acos, pow
from mathutils import Vector, Matrix, Euler
from bpy.types import Menu
from bpy.props import FloatProperty, BoolProperty, StringProperty, EnumProperty
from bpy_extras.object_utils import object_data_add
from add_mesh_extra_objects.add_mesh_solid import createPolys

bl_info = {
    'name': 'Dice Gen',
    'author': 'Long Tran',
    'version': (1, 0, 0),
    'blender': (2, 93, 0),
    'location': 'View3D > Add > Mesh',
    'description': 'Generate polyhedral dice models.',
    'category': 'Add Mesh',
    'doc_url': 'https://github.com/Longi94/blender-dice-gen/wiki',
    'tracker_url': 'https://github.com/Longi94/blender-dice-gen/issues'
}

NUMBER_IND_NONE = 'none'
NUMBER_IND_BAR = 'bar'
NUMBER_IND_PERIOD = 'period'

HALF_PI = math.pi / 2


def leg_b(leg_a, h):
    """given a leg of a right angle triangle and it's height, calculate the other leg"""
    return sqrt(pow(h, 2) + (pow(h, 4) / (pow(leg_a, 2) - pow(h, 2))))


# https://dmccooey.com/polyhedra
CONSTANTS = {
    'tetrahedron': {
        'dihedral_angle': acos(1 / 3),
        'height': sqrt(2 / 3),
        'c0': sqrt(2) / 4
    },
    'octahedron': {
        'dihedral_angle': acos(sqrt(5) / -5),
        'circumscribed_r': (sqrt(3) + sqrt(15)) / 4,
        'inscribed_r': sqrt(10 * (25 + 11 * sqrt(5))) / 20,
        'c0': (1 + sqrt(5)) / 4,
        'c1': (3 + sqrt(5)) / 4,
        'c2': 0.5
    },
    'icosahedron': {
        'dihedral_angle': acos(sqrt(5) / -3),
        'circumscribed_r': sqrt(2 * (5 + sqrt(5))) / 4,
        'inscribed_r': (3 * sqrt(3) + sqrt(15)) / 12,
        'c0': (1 + sqrt(5)) / 4,
        'c1': 0.5
    },
    'pentagonal_trap': {
        'inscribed_r': sqrt(5 * (5 + 2 * sqrt(5))) / 10,
        'base_height': 1.1180340051651,
        'base_width': leg_b(1.1180340051651, 0.5),
        'c0': (sqrt(5) - 1) / 4,
        'c1': (1 + sqrt(5)) / 4,
        'c2': (3 + sqrt(5)) / 4,
        'c3': 0.5
    }
}

# calculate rotation of trapezohedron to have it stand upright
# from dice-gen Math.acos((C0 - C2) / Math.sqrt(Math.pow(C0 - C2, 2) + 4 * Math.pow(C1, 2)))
CONSTANTS['pentagonal_trap']['angle'] = Euler((0, 0, acos(
    (CONSTANTS['pentagonal_trap']['c0'] - CONSTANTS['pentagonal_trap']['c2']) / sqrt(
        pow(CONSTANTS['pentagonal_trap']['c0'] - CONSTANTS['pentagonal_trap']['c2'], 2) + 4 * pow(
            CONSTANTS['pentagonal_trap']['c1'], 2)))), 'XYZ')

CONSTANTS['pentagonal_trap']['angle'].rotate(Euler((HALF_PI, 0, 0), 'XYZ'))


class Mesh:

    def __init__(self, name):
        self.vertices = None
        self.faces = None
        self.name = name
        self.dice_mesh = None
        self.base_font_scale = 1

    def create(self, context):
        self.dice_mesh = create_mesh(context, self.vertices, self.faces, self.name)
        # reset transforms
        self.dice_mesh.matrix_world = Matrix()
        return self.dice_mesh

    def get_numbers(self):
        return []

    def get_number_locations(self):
        return []

    def get_number_rotations(self):
        return []

    def create_numbers(self, context, size, number_scale, number_depth, font_path, one_offset,
                       number_indicator_type=NUMBER_IND_NONE, period_indicator_scale=1, period_indicator_space=1,
                       bar_indicator_height=1, bar_indicator_width=1, bar_indicator_space=1,
                       center_bar=True):
        numbers = self.get_numbers()
        locations = self.get_number_locations()
        rotations = self.get_number_rotations()

        font_size = self.base_font_scale * size * number_scale

        numbers_object = create_numbers(context, numbers, locations, rotations, font_path, font_size, number_depth,
                                        number_indicator_type, period_indicator_scale, period_indicator_space,
                                        bar_indicator_height, bar_indicator_width, bar_indicator_space,
                                        center_bar, one_offset)

        if numbers_object is not None:
            apply_boolean_modifier(self.dice_mesh, numbers_object)


class Tetrahedron(Mesh):
    def __init__(self, name, size, number_center_offset):
        super().__init__(name)
        self.size = size
        self.number_center_offset = number_center_offset

        c0 = CONSTANTS['tetrahedron']['c0'] / CONSTANTS['tetrahedron']['height'] * size

        self.vertices = [(c0, -c0, c0), (c0, c0, -c0), (-c0, c0, c0), (-c0, -c0, -c0)]
        self.faces = [[0, 1, 2], [1, 0, 3], [2, 3, 0], [3, 2, 1]]

        self.base_font_scale = 0.3

    def get_numbers(self):
        return [str(math.floor(i / 3) + 1) for i in range(12)]

    def get_number_locations(self):
        # face centers
        centers = [Vector((
            ((self.vertices[f[0]][0] + self.vertices[f[1]][0] + self.vertices[f[2]][0]) / 3),
            ((self.vertices[f[0]][1] + self.vertices[f[1]][1] + self.vertices[f[2]][1]) / 3),
            ((self.vertices[f[0]][2] + self.vertices[f[1]][2] + self.vertices[f[2]][2]) / 3)
        )) for f in self.faces]
        vertices = [Vector(v) for v in self.vertices]

        location_vectors = [
            centers[0].lerp(vertices[2], self.number_center_offset),
            centers[2].lerp(vertices[2], self.number_center_offset),
            centers[3].lerp(vertices[2], self.number_center_offset),
            centers[0].lerp(vertices[1], self.number_center_offset),
            centers[1].lerp(vertices[1], self.number_center_offset),
            centers[3].lerp(vertices[1], self.number_center_offset),
            centers[0].lerp(vertices[0], self.number_center_offset),
            centers[1].lerp(vertices[0], self.number_center_offset),
            centers[2].lerp(vertices[0], self.number_center_offset),
            centers[1].lerp(vertices[3], self.number_center_offset),
            centers[2].lerp(vertices[3], self.number_center_offset),
            centers[3].lerp(vertices[3], self.number_center_offset)
        ]

        return [(v.x, v.y, v.z) for v in location_vectors]

    def get_number_rotations(self):
        return [
            (CONSTANTS['tetrahedron']['dihedral_angle'] / 2, math.pi / 4, HALF_PI),
            (CONSTANTS['tetrahedron']['dihedral_angle'] / 2, -math.pi / 4, 0),
            ((math.pi - CONSTANTS['tetrahedron']['dihedral_angle']) / 2, 0, math.pi / 4),
            (-(math.pi - CONSTANTS['tetrahedron']['dihedral_angle']) / 2, 0, -math.pi / 4),
            (CONSTANTS['tetrahedron']['dihedral_angle'] / 2, math.pi * 3 / 4, 0),
            (CONSTANTS['tetrahedron']['dihedral_angle'] / 2, math.pi * 5 / 4, math.pi * 3 / 2),
            (CONSTANTS['tetrahedron']['dihedral_angle'] / 2, -math.pi / 4, -math.pi),
            (-(math.pi - CONSTANTS['tetrahedron']['dihedral_angle']) / 2, math.pi, -math.pi * 3 / 4),
            (CONSTANTS['tetrahedron']['dihedral_angle'] / 2, math.pi / 4, -HALF_PI),
            (CONSTANTS['tetrahedron']['dihedral_angle'] / 2, math.pi * 5 / 4, HALF_PI),
            (-(math.pi - CONSTANTS['tetrahedron']['dihedral_angle']) / 2, 0, math.pi * 3 / 4),
            (CONSTANTS['tetrahedron']['dihedral_angle'] / 2, math.pi * 3 / 4, math.pi)
        ]


class D4Crystal(Mesh):

    def __init__(self, name, size, base_height, point_height):
        super().__init__(name)
        self.size = size

        c0 = 0.5 * size
        c1 = 0.5 * base_height
        c2 = 0.5 * base_height + point_height

        self.vertices = [(-c0, -c0, c1), (c0, -c0, c1), (c0, c0, c1), (-c0, c0, c1), (-c0, -c0, -c1), (c0, -c0, -c1),
                         (c0, c0, -c1), (-c0, c0, -c1), (0, 0, c2), (0, 0, -c2)]
        self.faces = [[0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7], [0, 1, 8], [1, 2, 8], [2, 3, 8],
                      [3, 0, 8], [4, 5, 9], [5, 6, 9], [6, 7, 9], [7, 4, 9]]

        self.base_font_scale = 0.8

    def get_numbers(self):
        return numbers(4)

    def get_number_locations(self):
        c0 = 0.5 * self.size
        return [(c0, 0, 0), (0, c0, 0), (0, -c0, 0), (-c0, 0, 0)]

    def get_number_rotations(self):
        return [(HALF_PI, 0, HALF_PI), (HALF_PI, 0, HALF_PI * 2), (HALF_PI, 0, 0), (HALF_PI, 0, HALF_PI * 3)]


class D4Shard(Mesh):

    def __init__(self, name, size, top_point_height, bottom_point_height, number_v_offset):
        super().__init__(name)
        self.size = size
        self.number_v_offset = number_v_offset
        self.bottom_point_height = bottom_point_height

        c0 = size / sqrt(2)
        c1 = top_point_height * c0
        c2 = bottom_point_height * c0

        self.vertices = [(c0, 0, 0), (0, c0, 0), (0, -c0, 0), (-c0, 0, 0), (0, 0, c1), (0, 0, -c2)]
        self.faces = [[0, 1, 4], [1, 3, 4], [3, 2, 4], [2, 0, 4], [0, 1, 5], [1, 3, 5], [3, 2, 5], [2, 0, 5]]

        self.base_font_scale = 0.8

    def get_numbers(self):
        return numbers(4)

    def get_number_locations(self):
        c0 = self.size / 2 / sqrt(2) * self.number_v_offset
        c1 = self.size / sqrt(2) * self.bottom_point_height * (1 - self.number_v_offset)
        return [(c0, c0, -c1), (-c0, c0, -c1), (c0, -c0, -c1), (-c0, -c0, -c1)]

    def get_number_rotations(self):
        c0 = self.size / 2 / sqrt(2)
        c1 = self.size / sqrt(2) * self.bottom_point_height
        angle = math.pi / 2 + Vector((0, 0, c1)).angle(Vector((c0, c0, c1)))
        return [(angle, 0, math.pi * 3 / 4), (angle, 0, math.pi * 5 / 4), (angle, 0, math.pi * 1 / 4),
                (angle, 0, math.pi * 7 / 4)]


class Cube(Mesh):

    def __init__(self, name, size):
        super().__init__(name)

        # Calculate the necessary constants
        self.v_coord_const = 0.5 * size
        s = self.v_coord_const

        # create the vertices and faces
        self.vertices = [(-s, -s, -s), (s, -s, -s), (s, s, -s), (-s, s, -s), (-s, -s, s), (s, -s, s), (s, s, s),
                         (-s, s, s)]
        self.faces = [[0, 3, 2, 1], [0, 1, 5, 4], [0, 4, 7, 3], [6, 5, 1, 2], [6, 2, 3, 7], [6, 7, 4, 5]]

    def get_numbers(self):
        return numbers(6)

    def get_number_locations(self):
        s = self.v_coord_const
        return [(0, -s, 0), (-s, 0, 0), (0, 0, s), (0, 0, -s), (s, 0, 0), (0, s, 0)]

    def get_number_rotations(self):
        return [
            (HALF_PI, 0, 0),
            (math.pi, HALF_PI, 0),
            (0, 0, 0),
            (math.pi, 0, 0),
            (0, HALF_PI, 0),
            (-HALF_PI, 0, 0)
        ]


class Octahedron(Mesh):

    def __init__(self, name, size):
        super().__init__(name)

        # calculate circumscribed sphere radius from inscribed sphere radius
        # diameter of the inscribed sphere is the face 2 face length of the octahedron
        self.circumscribed_r = (size * math.sqrt(3)) / 2
        s = self.circumscribed_r

        # create the vertices and faces
        self.vertices = [(s, 0, 0), (-s, 0, 0), (0, s, 0), (0, -s, 0), (0, 0, s), (0, 0, -s)]
        self.faces = [[4, 0, 2], [4, 2, 1], [4, 1, 3], [4, 3, 0], [5, 2, 0], [5, 1, 2], [5, 3, 1], [5, 0, 3]]

        self.base_font_scale = 0.7

    def get_numbers(self):
        return numbers(8)

    def get_number_locations(self):
        c = self.circumscribed_r / 3
        return [
            (c, c, c),
            (c, c, -c),
            (c, -c, -c),
            (c, -c, c),
            (-c, c, -c),
            (-c, c, c),
            (-c, -c, c),
            (-c, -c, -c),
        ]

    def get_number_rotations(self):
        # dihedral angle / 2
        da = math.acos(-1 / 3)
        return [
            (da / 2, 0, math.pi * 3 / 4),
            (-math.pi + da / 2, 0, -math.pi / 4),
            (-math.pi + da / 2, 0, -math.pi * 3 / 4),
            (da / 2, 0, math.pi * 1 / 4),
            (-math.pi + da / 2, 0, math.pi * 1 / 4),
            (da / 2, 0, -math.pi * 3 / 4),
            (da / 2, 0, -math.pi * 1 / 4),
            (-math.pi + da / 2, 0, math.pi * 3 / 4),
        ]


class Dodecahedron(Mesh):

    def __init__(self, name, size):
        super().__init__(name)
        self.size = size

        # Calculate the necessary constants https://dmccooey.com/polyhedra/Dodecahedron.html
        edge_length = size / 2 / CONSTANTS['octahedron']['inscribed_r']

        c0 = CONSTANTS['octahedron']['c0'] * edge_length
        c1 = CONSTANTS['octahedron']['c1'] * edge_length
        s = CONSTANTS['octahedron']['c2'] * edge_length

        self.vertices = [(0.0, s, c1), (0.0, s, -c1), (0.0, -s, c1), (0.0, -s, -c1), (c1, 0.0, s), (c1, 0.0, -s),
                         (-c1, 0.0, s), (-c1, 0.0, -s), (s, c1, 0.0), (s, -c1, 0.0), (-s, c1, 0.0), (-s, -c1, 0.0),
                         (c0, c0, c0), (c0, c0, -c0), (c0, -c0, c0), (c0, -c0, -c0), (-c0, c0, c0), (-c0, c0, -c0),
                         (-c0, -c0, c0), (-c0, -c0, -c0)]

        self.faces = [[0, 2, 14, 4, 12], [0, 12, 8, 10, 16], [0, 16, 6, 18, 2], [7, 6, 16, 10, 17],
                      [7, 17, 1, 3, 19], [7, 19, 11, 18, 6], [9, 11, 19, 3, 15], [9, 15, 5, 4, 14],
                      [9, 14, 2, 18, 11], [13, 1, 17, 10, 8], [13, 8, 12, 4, 5], [13, 5, 15, 3, 1]]

        self.base_font_scale = 0.5

    def get_numbers(self):
        return numbers(12)

    def get_number_locations(self):
        dual_e = self.size / 2 / CONSTANTS['icosahedron']['circumscribed_r']
        c0 = dual_e * CONSTANTS['icosahedron']['c0']
        c1 = dual_e * CONSTANTS['icosahedron']['c1']
        return [
            (c1, 0.0, c0),
            (0.0, c0, c1),
            (-c1, 0.0, c0),
            (0.0, -c0, c1),
            (c0, -c1, 0.0),
            (c0, c1, 0.0),
            (-c0, -c1, 0.0),
            (-c0, c1, 0.0),
            (0.0, c0, -c1),
            (c1, 0.0, -c0),
            (0.0, -c0, -c1),
            (-c1, 0.0, -c0),
        ]

    def get_number_rotations(self):
        angles = [Euler((0, 0, 0), 'XYZ') for _ in range(12)]

        angles[0].z = math.radians(-162)
        angles[0].rotate(Euler((0, (math.pi - CONSTANTS['octahedron']['dihedral_angle']) / 2, 0), 'XYZ'))

        angles[1].z = math.radians(36)
        angles[1].rotate(Euler((CONSTANTS['octahedron']['dihedral_angle'] / -2, 0, 0), 'XYZ'))

        angles[2].z = HALF_PI
        angles[2].x = -(math.pi - CONSTANTS['octahedron']['dihedral_angle']) / 2

        angles[3].z = math.radians(144)
        angles[3].rotate(Euler((CONSTANTS['octahedron']['dihedral_angle'] / 2, 0, 0), 'XYZ'))

        angles[4].y = HALF_PI
        angles[4].rotate(Euler((-math.radians(108), 0, 0), 'XYZ'))
        angles[4].rotate(Euler((0, 0, (math.pi - CONSTANTS['octahedron']['dihedral_angle']) / -2), 'XYZ'))

        angles[5].y = HALF_PI
        angles[5].rotate(Euler((-math.radians(72), 0, 0), 'XYZ'))
        angles[5].rotate(Euler((0, 0, (math.pi - CONSTANTS['octahedron']['dihedral_angle']) / 2), 'XYZ'))

        angles[6].y = -HALF_PI
        angles[6].rotate(Euler((math.radians(108), 0, 0), 'XYZ'))
        angles[6].rotate(Euler((0, 0, (math.pi - CONSTANTS['octahedron']['dihedral_angle']) / 2), 'XYZ'))

        angles[7].y = HALF_PI
        angles[7].rotate(Euler((math.radians(72), 0, 0), 'XYZ'))
        angles[7].rotate(Euler((0, 0, -(math.pi - CONSTANTS['octahedron']['dihedral_angle']) / 2), 'XYZ'))

        angles[8].z = math.radians(-36)
        angles[8].y = math.pi
        angles[8].rotate(Euler((CONSTANTS['octahedron']['dihedral_angle'] / 2, 0, 0), 'XYZ'))

        angles[9].x = math.pi
        angles[9].z = HALF_PI
        angles[9].rotate(Euler((0, -(math.pi - CONSTANTS['octahedron']['dihedral_angle']) / 2, 0), 'XYZ'))

        angles[10].x = math.pi
        angles[10].z = math.radians(36)
        angles[10].rotate(Euler((-CONSTANTS['octahedron']['dihedral_angle'] / 2, 0, 0), 'XYZ'))

        angles[11].z = math.radians(342)
        angles[11].rotate(Euler((0, math.pi + (math.pi - CONSTANTS['octahedron']['dihedral_angle']) / 2, 0), 'XYZ'))

        return [(a.x, a.y, a.z) for a in angles]


class Icosahedron(Mesh):

    def __init__(self, name, size):
        super().__init__(name)
        self.size = size

        # Calculate the necessary constants https://dmccooey.com/polyhedra/Icosahedron.html
        edge_length = size / 2 / CONSTANTS['icosahedron']['inscribed_r']

        c0 = edge_length * CONSTANTS['icosahedron']['c0']
        c1 = edge_length * CONSTANTS['icosahedron']['c1']

        self.vertices = [(c1, 0.0, c0), (c1, 0.0, -c0), (-c1, 0.0, c0), (-c1, 0.0, -c0), (c0, c1, 0.0), (c0, -c1, 0.0),
                         (-c0, c1, 0.0), (-c0, -c1, 0.0), (0.0, c0, c1), (0.0, c0, -c1), (0.0, -c0, c1),
                         (0.0, -c0, -c1)]
        self.faces = [[0, 2, 10], [0, 10, 5], [0, 5, 4], [0, 4, 8], [0, 8, 2], [3, 1, 11], [3, 11, 7], [3, 7, 6],
                      [3, 6, 9], [3, 9, 1], [2, 6, 7], [2, 7, 10], [10, 7, 11], [10, 11, 5], [5, 11, 1], [5, 1, 4],
                      [4, 1, 9], [4, 9, 8], [8, 9, 6], [8, 6, 2]]

        self.base_font_scale = 0.3

    def get_numbers(self):
        return numbers(20)

    def get_number_locations(self):
        dual_e = self.size / 2 / CONSTANTS['octahedron']['circumscribed_r']

        c0 = CONSTANTS['octahedron']['c0'] * dual_e
        c1 = CONSTANTS['octahedron']['c1'] * dual_e
        s = CONSTANTS['octahedron']['c2'] * dual_e

        return [
            (0.0, s, c1),
            (-c0, -c0, -c0),
            (s, c1, 0.0),
            (s, -c1, 0.0),
            (-c0, -c0, c0),
            (c1, 0.0, -s),
            (-c0, c0, c0),
            (0.0, s, -c1),
            (c1, 0.0, s),
            (-c0, c0, -c0),
            (c0, -c0, c0),
            (-c1, 0.0, -s),
            (0.0, -s, c1),
            (c0, -c0, -c0),
            (-c1, 0.0, s),
            (c0, c0, -c0),
            (-s, c1, 0.0),
            (-s, -c1, 0.0),
            (c0, c0, c0),
            (0.0, -s, -c1)
        ]

    def get_number_rotations(self):
        angles = [Euler((0, 0, 0), 'XYZ') for _ in range(20)]

        # TODO magic numbers copied out of blender with alignment trick, try to find exact equation

        angles[0].x = -(math.pi - CONSTANTS['icosahedron']['dihedral_angle']) / 2

        angles[1].x = -0.918438
        angles[1].y = -2.82743
        angles[1].z = -4.15881

        angles[2].x = HALF_PI
        angles[2].y = 5 / 6 * math.pi
        angles[2].z = math.pi - ((math.pi - CONSTANTS['icosahedron']['dihedral_angle']) / 2)

        angles[3].x = HALF_PI
        angles[3].y = -math.pi / 6
        angles[3].z = (math.pi - CONSTANTS['icosahedron']['dihedral_angle']) / 2

        angles[4].x = -0.918438
        angles[4].y = -0.314159
        angles[4].z = 2.12437

        angles[5].x = HALF_PI
        angles[5].y = math.pi / 3
        angles[5].z = HALF_PI
        angles[5].rotate(Euler((0, (math.pi - CONSTANTS['icosahedron']['dihedral_angle']) / 2, 0), 'XYZ'))

        angles[6].x = -0.918438
        angles[6].y = 0.314159
        angles[6].z = 1.01722

        angles[7].x = -(math.pi - CONSTANTS['icosahedron']['dihedral_angle']) / 2
        angles[7].y = math.pi

        angles[8].x = -HALF_PI
        angles[8].y = -math.pi / 3
        angles[8].z = -HALF_PI
        angles[8].rotate(Euler((0, -(math.pi - CONSTANTS['icosahedron']['dihedral_angle']) / 2, 0), 'XYZ'))

        angles[9].x = -4.06003
        angles[9].y = 0.314159
        angles[9].z = -2.12437

        angles[10].x = -0.918438
        angles[10].y = 0.314159
        angles[10].z = -2.12437

        angles[11].x = HALF_PI
        angles[11].y = -math.pi / 3
        angles[11].z = -HALF_PI
        angles[11].rotate(Euler((0, -(math.pi - CONSTANTS['icosahedron']['dihedral_angle']) / 2, 0), 'XYZ'))

        angles[12].x = -(math.pi - CONSTANTS['icosahedron']['dihedral_angle']) / 2
        angles[12].z = math.pi

        angles[13].x = 2.22315
        angles[13].y = 0.314159
        angles[13].z = 1.01722

        angles[14].x = -HALF_PI
        angles[14].y = math.pi / 3
        angles[14].z = HALF_PI
        angles[14].rotate(Euler((0, (math.pi - CONSTANTS['icosahedron']['dihedral_angle']) / 2, 0), 'XYZ'))

        angles[15].x = -0.918438
        angles[15].y = -2.82743
        angles[15].z = -1.01722

        angles[16].x = HALF_PI
        angles[16].y = 7 / 6 * math.pi
        angles[16].z = math.pi + ((math.pi - CONSTANTS['icosahedron']['dihedral_angle']) / 2)

        angles[17].x = HALF_PI
        angles[17].y = math.pi / 6
        angles[17].z = -(math.pi - CONSTANTS['icosahedron']['dihedral_angle']) / 2

        angles[18].x = -0.918438
        angles[18].y = -0.314159
        angles[18].z = -1.01722

        angles[19].x = math.pi
        angles[19].z = 2 / 3 * math.pi
        angles[19].rotate(Euler((-(math.pi - CONSTANTS['icosahedron']['dihedral_angle']) / 2, 0, 0), 'XYZ'))

        return [(a.x, a.y, a.z) for a in angles]


class SquashedPentagonalTrapezohedron(Mesh):

    def __init__(self, name, size, height, number_v_offset):
        super().__init__(name)
        self.size = size
        self.height = height
        self.number_v_offset = number_v_offset

        antiprism_e = size / 2 / CONSTANTS['pentagonal_trap']['inscribed_r']

        c0 = CONSTANTS['pentagonal_trap']['c0'] * antiprism_e
        c1 = CONSTANTS['pentagonal_trap']['c1'] * antiprism_e
        c2 = CONSTANTS['pentagonal_trap']['c2'] * antiprism_e
        c3 = CONSTANTS['pentagonal_trap']['c3'] * antiprism_e

        scaled_base_height = CONSTANTS['pentagonal_trap']['base_height'] * size
        scaled_base_width = CONSTANTS['pentagonal_trap']['base_width'] * size

        scaled_height = scaled_base_height * height
        scaled_width = leg_b(scaled_height, size / 2)
        width = scaled_width / scaled_base_width

        # TODO figure out where this angle comes from
        self.vertices = [(0.0, c0, c1), (0.0, c0, -c1), (0.0, -c0, c1), (0.0, -c0, -c1), (c3, c3, c3), (c3, c3, -c3),
                         (-c3, -c3, c3), (-c3, -c3, -c3), (c2, -c1, 0.0), (-c2, c1, 0.0), (c0, c1, 0.0),
                         (-c0, -c1, 0.0)]

        def transform(v):
            # rotate the vectors, so the trapezohedron is up right
            vector = Vector(v)
            vector.rotate(CONSTANTS['pentagonal_trap']['angle'])

            # scale the body
            vector.z *= height
            vector.y *= width
            vector.x *= width
            return vector.x, vector.y, vector.z

        self.vertices = list(map(transform, self.vertices))

        self.faces = [[8, 2, 6, 11], [8, 11, 7, 3], [8, 3, 1, 5], [8, 5, 10, 4], [8, 4, 0, 2], [9, 0, 4, 10],
                      [9, 10, 5, 1], [9, 1, 3, 7], [9, 7, 11, 6], [9, 6, 2, 0]]

    def get_number_locations(self):
        vectors = [Vector(v) for v in self.vertices]

        lerp_factor = self.number_v_offset
        location_vectors = [
            vectors[6].lerp(vectors[8], lerp_factor),
            vectors[3].lerp(vectors[9], lerp_factor),
            vectors[1].lerp(vectors[8], lerp_factor),
            vectors[4].lerp(vectors[9], lerp_factor),
            vectors[10].lerp(vectors[8], lerp_factor),
            vectors[11].lerp(vectors[9], lerp_factor),
            vectors[7].lerp(vectors[8], lerp_factor),
            vectors[2].lerp(vectors[9], lerp_factor),
            vectors[0].lerp(vectors[8], lerp_factor),
            vectors[5].lerp(vectors[9], lerp_factor)
        ]

        return [(v.x, v.y, v.z) for v in location_vectors]

    def get_number_rotations(self):
        a = Vector(self.vertices[9])
        b = Vector(self.vertices[10]) - Vector(self.vertices[8])
        number_angle = HALF_PI - a.angle(b)
        return [
            (number_angle, 0, -HALF_PI - math.pi * 6 / 5),
            (math.pi + number_angle, 0, -HALF_PI - math.pi * 8 / 5),
            (number_angle, 0, -HALF_PI - math.pi * 2 / 5),
            (math.pi + number_angle, 0, -HALF_PI - math.pi * 4 / 5),
            (number_angle, 0, -HALF_PI),
            (math.pi + number_angle, 0, -HALF_PI),
            (number_angle, 0, -HALF_PI - math.pi * 4 / 5),
            (math.pi + number_angle, 0, -HALF_PI - math.pi * 2 / 5),
            (number_angle, 0, -HALF_PI - math.pi * 8 / 5),
            (math.pi + number_angle, 0, -HALF_PI - math.pi * 6 / 5)
        ]


class D10Mesh(SquashedPentagonalTrapezohedron):

    def __init__(self, name, size, height, number_v_offset):
        super().__init__(name, size, height, number_v_offset)
        self.base_font_scale = 0.6

    def get_numbers(self):
        return [str((i + 1) % 10) for i in range(10)]


class D100Mesh(SquashedPentagonalTrapezohedron):

    def __init__(self, name, size, height, number_v_offset):
        super().__init__(name, size, height, number_v_offset)
        self.base_font_scale = 0.45

    def get_numbers(self):
        return [f'{str((i + 1) % 10)}0' for i in range(10)]


def numbers(n: int) -> List[str]:
    return [str(i + 1) for i in range(n)]


def set_origin(o, v):
    """
    set origin to a specific location
    """
    me = o.data
    mw = o.matrix_world
    current = o.location
    T = Matrix.Translation(current - v)
    me.transform(T)
    mw.translation = mw @ v


def set_origin_center_bounds(o):
    """
    set an objects origin to the center of its bounding box
    :param o:
    :return:
    """
    me = o.data

    max_x = max((v.co.x for v in me.vertices))
    max_y = max((v.co.y for v in me.vertices))
    max_z = max((v.co.z for v in me.vertices))

    min_x = min((v.co.x for v in me.vertices))
    min_y = min((v.co.y for v in me.vertices))
    min_z = min((v.co.z for v in me.vertices))

    set_origin(o, Vector(((min_x + max_x) / 2, (min_y + max_y) / 2, (min_z + max_z) / 2)))


def set_origin_min_bounds(o):
    """
    set an objects origin to the bottom_left corner of its bounding box
    :param o:
    :return:
    """
    me = o.data

    max_z = max((v.co.z for v in me.vertices))
    min_x = min((v.co.x for v in me.vertices))
    min_y = min((v.co.y for v in me.vertices))
    min_z = min((v.co.z for v in me.vertices))

    set_origin(o, Vector((min_x, min_y, (min_z + max_z) / 2)))


def create_mesh(context, vertices, faces, name):
    verts = [Vector(i) for i in vertices]

    # turn n-gons in quads and tri's
    faces = createPolys(faces)

    # generate object
    # Create new mesh
    mesh = bpy.data.meshes.new(name)

    # Make a mesh from a list of verts/edges/faces.
    mesh.from_pydata(verts, [], faces)

    # Update mesh geometry after adding stuff.
    mesh.update()

    return object_data_add(context, mesh, operator=None)


def apply_transform(ob, use_location=False, use_rotation=False, use_scale=False):
    """
    https://blender.stackexchange.com/questions/159538/how-to-apply-all-transformations-to-an-object-at-low-level
    :param ob:
    :param use_location:
    :param use_rotation:
    :param use_scale:
    :return:
    """
    mb = ob.matrix_basis
    I = Matrix()
    loc, rot, scale = mb.decompose()

    # rotation
    T = Matrix.Translation(loc)
    # R = rot.to_matrix().to_4x4()
    R = mb.to_3x3().normalized().to_4x4()
    S = Matrix.Diagonal(scale).to_4x4()

    transform = [I, I, I]
    basis = [T, R, S]

    def swap(i):
        transform[i], basis[i] = basis[i], transform[i]

    if use_location:
        swap(0)
    if use_rotation:
        swap(1)
    if use_scale:
        swap(2)

    M = transform[0] @ transform[1] @ transform[2]
    if hasattr(ob.data, "transform"):
        ob.data.transform(M)
    for c in ob.children:
        c.matrix_local = M @ c.matrix_local

    ob.matrix_basis = basis[0] @ basis[1] @ basis[2]


def join(objects):
    bpy.context.view_layer.objects.active = objects[0]
    ctx = bpy.context.copy()
    # one of the objects to join
    ctx['active_object'] = objects[0]
    ctx['selected_editable_objects'] = objects
    # TODO apparently, ops calls trigger a viewport refresh, should find a way to replace with a low level join
    bpy.ops.object.join(ctx)
    return objects[0]


def validate_font_path(filepath):
    # set font to emtpy if it's not a ttf file
    if filepath and os.path.splitext(filepath)[1].lower() not in ('.ttf', '.otf'):
        return ''
    return filepath


def get_font(filepath):
    if filepath:
        bpy.data.fonts.load(filepath=filepath, check_existing=True)
        return next(filter(lambda x: x.filepath == filepath, bpy.data.fonts))
    else:
        bpy.data.fonts.load(filepath='<builtin>', check_existing=True)
        return bpy.data.fonts['Bfont']


def apply_boolean_modifier(body_object, numbers_object):
    """
    Add a BOOLEAN modifier to body_object that targets
    :param context:
    :param body_object:
    :param numbers_object
    :return:
    """
    numbers_boolean = body_object.modifiers.new(type='BOOLEAN', name='boolean')
    numbers_boolean.object = bpy.data.objects[numbers_object.name]
    numbers_boolean.show_viewport = False


def create_text_mesh(context, text, font_path, font_size, name, extrude=0):
    # load the font
    font = get_font(font_path)

    # create the text curve
    font_curve = bpy.data.curves.new(type='FONT', name=name)
    font_curve.body = text
    font_curve.font = font
    font_curve.size = font_size
    font_curve.extrude = extrude
    font_curve.offset = 0

    # create object from curve
    curve_obj = bpy.data.objects.new('temp_curve_obj', font_curve)

    # convert curve to mesh
    mesh = curve_obj.to_mesh().copy()
    curve_obj.to_mesh_clear()
    bpy.data.objects.remove(curve_obj)
    bpy.data.curves.remove(font_curve)
    return object_data_add(context, mesh, operator=None)


def create_numbers(context, numbers, locations, rotations, font_path, font_size, number_depth, number_indicator_type,
                   period_indicator_scale, period_indicator_space, bar_indicator_height, bar_indicator_width,
                   bar_indicator_space, center_bar, one_offset):
    number_objs = []
    # create the number meshes
    for i in range(len(locations)):
        number_object = create_number(context, numbers[i], font_path, font_size, number_depth, locations[i],
                                      rotations[i], number_indicator_type, period_indicator_scale,
                                      period_indicator_space, bar_indicator_height, bar_indicator_width,
                                      bar_indicator_space, center_bar, one_offset)
        number_objs.append(number_object)

    # join the numbers into a single object
    if len(number_objs):
        numbers = join(number_objs)
        apply_transform(numbers, use_rotation=True, use_location=True)
        return numbers

    return None


def create_number(context, number, font_path, font_size, number_depth, location, rotation, number_indicator_type,
                  period_indicator_scale, period_indicator_space, bar_indicator_height, bar_indicator_width,
                  bar_indicator_space, center_bar, one_offset):
    """
    Create a number mesh that will be used in a boolean modifier
    """
    # add number
    mesh_object = create_text_mesh(context, number, font_path, font_size, f'number_{number}', number_depth)

    # set origin to bounding box center
    set_origin_center_bounds(mesh_object)

    if number == '1':
        if one_offset > 0:
            number_width = mesh_object.dimensions.x
            new_origin = Vector((mesh_object.location.x + number_width * one_offset, mesh_object.location.y,
                                 mesh_object.location.z))
            set_origin(mesh_object, new_origin)
            pass
    elif number in ('6', '9'):
        if number_indicator_type == NUMBER_IND_PERIOD:
            p_obj = create_text_mesh(context, '.', font_path, font_size * period_indicator_scale, f'period_{number}',
                                     number_depth)

            # move origin of period to the bottom left corner of the mesh
            set_origin_min_bounds(p_obj)

            space = (1 / 20) * font_size * period_indicator_space

            # move period to the bottom right of the number
            p_obj.location = Vector((mesh_object.location.x + (mesh_object.dimensions.x / 2) + space,
                                     mesh_object.location.y - (mesh_object.dimensions.y / 2), 0))

            # join the period to the number
            mesh_object = join([mesh_object, p_obj])
        elif number_indicator_type == NUMBER_IND_BAR:
            # create a simple rectangle
            bar_width = mesh_object.dimensions.x * bar_indicator_width
            bar_height = (1 / 15) * font_size * bar_indicator_height
            bar_space = (1 / 20) * font_size * bar_indicator_space
            bar_obj = create_mesh(context,
                                  [(-bar_width / 2, -bar_space, number_depth),
                                   (bar_width / 2, -bar_space, number_depth),
                                   (-bar_width / 2, -bar_space - bar_height, number_depth),
                                   (bar_width / 2, -bar_space - bar_height, number_depth),
                                   (-bar_width / 2, -bar_space, -number_depth),
                                   (bar_width / 2, -bar_space, -number_depth),
                                   (-bar_width / 2, -bar_space - bar_height, -number_depth),
                                   (bar_width / 2, -bar_space - bar_height, -number_depth)],
                                  [[0, 1, 3, 2], [2, 3, 7, 6], [3, 1, 5, 7], [1, 0, 4, 5], [0, 2, 6, 4], [4, 6, 7, 5]],
                                  'bar_indicator')

            # move bar below the number
            bar_obj.location = Vector(
                (mesh_object.location.x, mesh_object.location.y - (mesh_object.dimensions.y / 2), 0))

            # join the bar to the number
            mesh_object = join([mesh_object, bar_obj])

            # recenter the mesh
            if center_bar:
                mesh_object.location = Vector((0, 0, 0))
                set_origin_center_bounds(mesh_object)

    mesh_object.location.x = location[0]
    mesh_object.location.y = location[1]
    mesh_object.location.z = location[2]

    mesh_object.rotation_euler.x = rotation[0]
    mesh_object.rotation_euler.y = rotation[1]
    mesh_object.rotation_euler.z = rotation[2]

    for f in mesh_object.data.polygons:
        f.use_smooth = False

    return mesh_object


def execute_generator(op, context, mesh_cls, name, **kwargs):
    # set font to emtpy if it's not a ttf file
    op.font_path = validate_font_path(op.font_path)

    # create the cube mesh
    die = mesh_cls(name, op.size, **kwargs)
    die.create(context)

    # create number curves
    if op.add_numbers:
        if op.number_indicator_type == NUMBER_IND_NONE:
            die.create_numbers(context, op.size, op.number_scale, op.number_depth, op.font_path, op.one_offset)
        else:
            die.create_numbers(context, op.size, op.number_scale, op.number_depth, op.font_path, op.one_offset,
                               op.number_indicator_type, op.period_indicator_scale, op.period_indicator_space,
                               op.bar_indicator_height, op.bar_indicator_width, op.bar_indicator_space, op.center_bar)

    return {'FINISHED'}


class D4Generator(bpy.types.Operator):
    """Generate a D4"""
    bl_idname = 'mesh.d4_add'
    bl_label = 'D4 Tetrahedron'
    bl_description = 'Generate a tetrahedron dice'
    bl_options = {'REGISTER', 'UNDO'}

    number_indicator_type = NUMBER_IND_NONE

    size: FloatProperty(
        name='Face2Face Length',
        description='Face-to-face size of the die',
        min=1,
        soft_min=1,
        max=100,
        soft_max=100,
        default=20,
        unit='LENGTH'
    )

    add_numbers: BoolProperty(
        name='Generate Numbers',
        default=True
    )

    number_scale: FloatProperty(
        name='Number Scale',
        description='Size of the numbers on the die',
        min=0.1,
        soft_min=0.1,
        max=2,
        soft_max=2,
        default=1
    )

    number_depth: FloatProperty(
        name='Number Depth',
        description='Depth of the numbers on the die',
        min=0.1,
        soft_min=0.1,
        max=2,
        soft_max=2,
        default=0.75,
        unit='LENGTH'
    )

    font_path: StringProperty(
        name='Font',
        description='Number font',
        maxlen=1024,
        subtype='FILE_PATH'
    )

    number_center_offset: FloatProperty(
        name='Number Center Offset',
        description='Distance of numbers from the center of a face',
        min=0.0,
        soft_min=0.0,
        max=1,
        soft_max=1,
        default=0.5
    )

    one_offset: FloatProperty(
        name='Number 1 Offset',
        description='Offset the number 1 horizontally for an alternative centering',
        min=0,
        soft_min=0,
        max=1,
        soft_max=1,
        default=0
    )

    def execute(self, context):
        return execute_generator(self, context, Tetrahedron, 'd4', number_center_offset=self.number_center_offset)


class D4CrystalGenerator(bpy.types.Operator):
    """Generate a D4 crystal"""
    bl_idname = 'mesh.d4_crystal_add'
    bl_label = 'D4 Crystal'
    bl_description = 'Generate a D4 crystal dice'
    bl_options = {'REGISTER', 'UNDO'}

    number_indicator_type = NUMBER_IND_NONE

    size: FloatProperty(
        name='Face2point Length',
        description='Face-to-point size of the die',
        min=1,
        soft_min=1,
        max=100,
        soft_max=100,
        default=12,
        unit='LENGTH'
    )

    base_height: FloatProperty(
        name='Base Height',
        description='Base height of the die (height of a face)',
        min=1,
        soft_min=1,
        max=100,
        soft_max=100,
        default=14,
        unit='LENGTH'
    )

    point_height: FloatProperty(
        name='Point Height',
        description='Point height of the die',
        min=1,
        soft_min=1,
        max=100,
        soft_max=100,
        default=7,
        unit='LENGTH'
    )

    add_numbers: BoolProperty(
        name='Generate Numbers',
        default=True
    )

    number_scale: FloatProperty(
        name='Number Scale',
        description='Size of the numbers on the die',
        min=0.1,
        soft_min=0.1,
        max=2,
        soft_max=2,
        default=1
    )

    number_depth: FloatProperty(
        name='Number Depth',
        description='Depth of the numbers on the die',
        min=0.1,
        soft_min=0.1,
        max=2,
        soft_max=2,
        default=0.75,
        unit='LENGTH'
    )

    font_path: StringProperty(
        name='Font',
        description='Number font',
        maxlen=1024,
        subtype='FILE_PATH'
    )

    one_offset: FloatProperty(
        name='Number 1 Offset',
        description='Offset the number 1 horizontally for an alternative centering',
        min=0,
        soft_min=0,
        max=1,
        soft_max=1,
        default=0
    )

    def execute(self, context):
        return execute_generator(self, context, D4Crystal, 'd4Crystal', base_height=self.base_height,
                                 point_height=self.point_height)


class D4ShardGenerator(bpy.types.Operator):
    """Generate a D4 shard"""
    bl_idname = 'mesh.d4_shard_add'
    bl_label = 'D4 Shard'
    bl_description = 'Generate a D4 crystal dice'
    bl_options = {'REGISTER', 'UNDO'}

    number_indicator_type = NUMBER_IND_NONE

    size: FloatProperty(
        name='Edge2edge length',
        description='Distance between 2 opposite horizontal edges',
        min=1,
        soft_min=1,
        max=100,
        soft_max=100,
        default=12,
        unit='LENGTH'
    )

    top_point_height: FloatProperty(
        name='Top Point Height',
        description='Top point height of the die',
        min=0.25,
        soft_min=0.25,
        max=2,
        soft_max=2,
        default=0.75
    )

    bottom_point_height: FloatProperty(
        name='Bottom Point Height',
        description='Bottom point height of the die',
        min=0.25,
        soft_min=0.25,
        max=2.5,
        soft_max=2.5,
        default=1.75
    )

    add_numbers: BoolProperty(
        name='Generate Numbers',
        default=True
    )

    number_scale: FloatProperty(
        name='Number Scale',
        description='Size of the numbers on the die',
        min=0.1,
        soft_min=0.1,
        max=2,
        soft_max=2,
        default=1
    )

    number_depth: FloatProperty(
        name='Number Depth',
        description='Depth of the numbers on the die',
        min=0.1,
        soft_min=0.1,
        max=2,
        soft_max=2,
        default=0.75,
        unit='LENGTH'
    )

    font_path: StringProperty(
        name='Font',
        description='Number font',
        maxlen=1024,
        subtype='FILE_PATH'
    )

    number_v_offset: FloatProperty(
        name='Number V Offset',
        description='Vertical offset of the number positioning',
        min=0.0,
        soft_min=0.0,
        max=1,
        soft_max=1,
        default=0.7
    )

    one_offset: FloatProperty(
        name='Number 1 Offset',
        description='Offset the number 1 horizontally for an alternative centering',
        min=0,
        soft_min=0,
        max=1,
        soft_max=1,
        default=0
    )

    def execute(self, context):
        return execute_generator(self, context, D4Shard, 'd4Shard', top_point_height=self.top_point_height,
                                 bottom_point_height=self.bottom_point_height, number_v_offset=self.number_v_offset)


class D6Generator(bpy.types.Operator):
    """Generate a D6"""
    bl_idname = 'mesh.d6_add'
    bl_label = 'D6 Cube'
    bl_description = 'Generate a cube dice'
    bl_options = {'REGISTER', 'UNDO'}

    size: FloatProperty(
        name='Face2face Length',
        description='Face-to-face size of the die',
        min=1,
        soft_min=1,
        max=100,
        soft_max=100,
        default=16,
        unit='LENGTH'
    )

    add_numbers: BoolProperty(
        name='Generate Numbers',
        default=True
    )

    number_scale: FloatProperty(
        name='Number Scale',
        description='Size of the numbers on the die',
        min=0.1,
        soft_min=0.1,
        max=2,
        soft_max=2,
        default=1
    )

    number_depth: FloatProperty(
        name='Number Depth',
        description='Depth of the numbers on the die',
        min=0.1,
        soft_min=0.1,
        max=2,
        soft_max=2,
        default=0.75,
        unit='LENGTH'
    )

    font_path: StringProperty(
        name='Font',
        description='Number font',
        maxlen=1024,
        subtype='FILE_PATH'
    )

    number_indicator_type: EnumProperty(
        name='Orientation Indicator',
        items=(
            (NUMBER_IND_NONE, 'None', ','),
            (NUMBER_IND_BAR, 'Bar', ''),
            (NUMBER_IND_PERIOD, 'Period', ''),
        ),
        default=NUMBER_IND_NONE,
        description='Orientation indicator for numbers 6 and 9'
    )

    period_indicator_scale: FloatProperty(
        name='Period Scale',
        description='Scale of the period orientation indicator',
        min=0.1,
        soft_min=0.1,
        max=2,
        soft_max=2,
        default=1
    )

    period_indicator_space: FloatProperty(
        name='Period Space',
        description='Space between the period orientation indicator and the number',
        min=0,
        soft_min=0,
        max=3,
        soft_max=3,
        default=1
    )

    bar_indicator_height: FloatProperty(
        name='Bar Height',
        description='Height scale of the bar orientation indicator',
        min=0.1,
        soft_min=0.1,
        max=3,
        soft_max=3,
        default=1
    )

    bar_indicator_width: FloatProperty(
        name='Bar Width',
        description='Width scale of the bar orientation indicator',
        min=0.1,
        soft_min=0.1,
        max=2,
        soft_max=2,
        default=1
    )

    bar_indicator_space: FloatProperty(
        name='Bar Space',
        description='Space between the bar orientation indicator and the number',
        min=0,
        soft_min=0,
        max=3,
        soft_max=3,
        default=1
    )

    center_bar: BoolProperty(
        name='Center Align Bar',
        description='If true, the bar indicator is included in the vertical alignment of the number',
        default=True
    )

    one_offset: FloatProperty(
        name='Number 1 Offset',
        description='Offset the number 1 horizontally for an alternative centering',
        min=0,
        soft_min=0,
        max=1,
        soft_max=1,
        default=0
    )

    def execute(self, context):
        return execute_generator(self, context, Cube, 'd6')


class D8Generator(bpy.types.Operator):
    """Generate a D8"""
    bl_idname = 'mesh.d8_add'
    bl_label = 'D8 Octahedron'
    bl_description = 'Generate a octahedron dice'
    bl_options = {'REGISTER', 'UNDO'}

    size: FloatProperty(
        name='Face2face Length',
        description='Face-to-face size of the die',
        min=1,
        soft_min=1,
        max=100,
        soft_max=100,
        default=15,
        unit='LENGTH'
    )

    add_numbers: BoolProperty(
        name='Generate Numbers',
        default=True
    )

    number_scale: FloatProperty(
        name='Number Scale',
        description='Size of the numbers on the die',
        min=0.1,
        soft_min=0.1,
        max=2,
        soft_max=2,
        default=1
    )

    number_depth: FloatProperty(
        name='Number Depth',
        description='Depth of the numbers on the die',
        min=0.1,
        soft_min=0.1,
        max=2,
        soft_max=2,
        default=0.75,
        unit='LENGTH'
    )

    font_path: StringProperty(
        name='Font',
        description='Number font',
        maxlen=1024,
        subtype='FILE_PATH'
    )

    number_indicator_type: EnumProperty(
        name='Orientation Indicator',
        items=(
            (NUMBER_IND_NONE, 'None', ','),
            (NUMBER_IND_BAR, 'Bar', ''),
            (NUMBER_IND_PERIOD, 'Period', ''),
        ),
        default=NUMBER_IND_NONE,
        description='Orientation indicator for numbers 6 and 9'
    )

    period_indicator_scale: FloatProperty(
        name='Period Scale',
        description='Scale of the period orientation indicator',
        min=0.1,
        soft_min=0.1,
        max=2,
        soft_max=2,
        default=1
    )

    period_indicator_space: FloatProperty(
        name='Period Space',
        description='Space between the period orientation indicator and the number',
        min=0,
        soft_min=0,
        max=3,
        soft_max=3,
        default=1
    )

    bar_indicator_height: FloatProperty(
        name='Bar Height',
        description='Height scale of the bar orientation indicator',
        min=0.1,
        soft_min=0.1,
        max=3,
        soft_max=3,
        default=1
    )

    bar_indicator_width: FloatProperty(
        name='Bar Width',
        description='Width scale of the bar orientation indicator',
        min=0.1,
        soft_min=0.1,
        max=2,
        soft_max=2,
        default=1
    )

    bar_indicator_space: FloatProperty(
        name='Bar Space',
        description='Space between the bar orientation indicator and the number',
        min=0,
        soft_min=0,
        max=3,
        soft_max=3,
        default=1
    )

    center_bar: BoolProperty(
        name='Center Align Bar',
        description='If true, the bar indicator is included in the vertical alignment of the number',
        default=True
    )

    one_offset: FloatProperty(
        name='Number 1 Offset',
        description='Offset the number 1 horizontally for an alternative centering',
        min=0,
        soft_min=0,
        max=1,
        soft_max=1,
        default=0
    )

    def execute(self, context):
        return execute_generator(self, context, Octahedron, 'd8')


class D12Generator(bpy.types.Operator):
    """Generate a D12"""
    bl_idname = 'mesh.d12_add'
    bl_label = 'D12 Dodecahedron'
    bl_description = 'Generate a dodecahedron dice'
    bl_options = {'REGISTER', 'UNDO'}

    size: FloatProperty(
        name='Face2face Length',
        description='Face-to-face size of the die',
        min=1,
        soft_min=1,
        max=100,
        soft_max=100,
        default=18,
        unit='LENGTH'
    )

    add_numbers: BoolProperty(
        name='Generate Numbers',
        default=True
    )

    number_scale: FloatProperty(
        name='Number Scale',
        description='Size of the numbers on the die',
        min=0.1,
        soft_min=0.1,
        max=2,
        soft_max=2,
        default=1
    )

    number_depth: FloatProperty(
        name='Number Depth',
        description='Depth of the numbers on the die',
        min=0.1,
        soft_min=0.1,
        max=2,
        soft_max=2,
        default=0.75,
        unit='LENGTH'
    )

    font_path: StringProperty(
        name='Font',
        description='Number font',
        maxlen=1024,
        subtype='FILE_PATH'
    )

    number_indicator_type: EnumProperty(
        name='Orientation Indicator',
        items=(
            (NUMBER_IND_NONE, 'None', ','),
            (NUMBER_IND_BAR, 'Bar', ''),
            (NUMBER_IND_PERIOD, 'Period', ''),
        ),
        default=NUMBER_IND_PERIOD,
        description='Orientation indicator for numbers 6 and 9'
    )

    period_indicator_scale: FloatProperty(
        name='Period Scale',
        description='Scale of the period orientation indicator',
        min=0.1,
        soft_min=0.1,
        max=2,
        soft_max=2,
        default=1
    )

    period_indicator_space: FloatProperty(
        name='Period Space',
        description='Space between the period orientation indicator and the number',
        min=0,
        soft_min=0,
        max=3,
        soft_max=3,
        default=1
    )

    bar_indicator_height: FloatProperty(
        name='Bar Height',
        description='Height scale of the bar orientation indicator',
        min=0.1,
        soft_min=0.1,
        max=3,
        soft_max=3,
        default=1
    )

    bar_indicator_width: FloatProperty(
        name='Bar Width',
        description='Width scale of the bar orientation indicator',
        min=0.1,
        soft_min=0.1,
        max=2,
        soft_max=2,
        default=1
    )

    bar_indicator_space: FloatProperty(
        name='Bar Space',
        description='Space between the bar orientation indicator and the number',
        min=0,
        soft_min=0,
        max=3,
        soft_max=3,
        default=1
    )

    center_bar: BoolProperty(
        name='Center Align Bar',
        description='If true, the bar indicator is included in the vertical alignment of the number',
        default=True
    )

    one_offset: FloatProperty(
        name='Number 1 Offset',
        description='Offset the number 1 horizontally for an alternative centering',
        min=0,
        soft_min=0,
        max=1,
        soft_max=1,
        default=0
    )

    def execute(self, context):
        return execute_generator(self, context, Dodecahedron, 'd12')


class D20Generator(bpy.types.Operator):
    """Generate a D20"""
    bl_idname = 'mesh.d20_add'
    bl_label = 'D20 Icosahedron'
    bl_description = 'Generate an icosahedron dice'
    bl_options = {'REGISTER', 'UNDO'}

    size: FloatProperty(
        name='Face2face Length',
        description='Face-to-face size of the die',
        min=1,
        soft_min=1,
        max=100,
        soft_max=100,
        default=20,
        unit='LENGTH'
    )

    add_numbers: BoolProperty(
        name='Generate Numbers',
        default=True
    )

    number_scale: FloatProperty(
        name='Number Scale',
        description='Size of the numbers on the die',
        min=0.1,
        soft_min=0.1,
        max=2,
        soft_max=2,
        default=1
    )

    number_depth: FloatProperty(
        name='Number Depth',
        description='Depth of the numbers on the die',
        min=0.1,
        soft_min=0.1,
        max=2,
        soft_max=2,
        default=0.75,
        unit='LENGTH'
    )

    font_path: StringProperty(
        name='Font',
        description='Number font',
        maxlen=1024,
        subtype='FILE_PATH'
    )

    number_indicator_type: EnumProperty(
        name='Orientation Indicator',
        items=(
            (NUMBER_IND_NONE, 'None', ','),
            (NUMBER_IND_BAR, 'Bar', ''),
            (NUMBER_IND_PERIOD, 'Period', ''),
        ),
        default=NUMBER_IND_PERIOD,
        description='Orientation indicator for numbers 6 and 9'
    )

    period_indicator_scale: FloatProperty(
        name='Period Scale',
        description='Scale of the period orientation indicator',
        min=0.1,
        soft_min=0.1,
        max=2,
        soft_max=2,
        default=1
    )

    period_indicator_space: FloatProperty(
        name='Period Space',
        description='Space between the period orientation indicator and the number',
        min=0,
        soft_min=0,
        max=3,
        soft_max=3,
        default=1
    )

    bar_indicator_height: FloatProperty(
        name='Bar Height',
        description='Height scale of the bar orientation indicator',
        min=0.1,
        soft_min=0.1,
        max=3,
        soft_max=3,
        default=1
    )

    bar_indicator_width: FloatProperty(
        name='Bar Width',
        description='Width scale of the bar orientation indicator',
        min=0.1,
        soft_min=0.1,
        max=2,
        soft_max=2,
        default=1
    )

    bar_indicator_space: FloatProperty(
        name='Bar Space',
        description='Space between the bar orientation indicator and the number',
        min=0,
        soft_min=0,
        max=3,
        soft_max=3,
        default=1
    )

    center_bar: BoolProperty(
        name='Center Align Bar',
        description='If true, the bar indicator is included in the vertical alignment of the number',
        default=True
    )

    one_offset: FloatProperty(
        name='Number 1 Offset',
        description='Offset the number 1 horizontally for an alternative centering',
        min=0,
        soft_min=0,
        max=1,
        soft_max=1,
        default=0
    )

    def execute(self, context):
        return execute_generator(self, context, Icosahedron, 'd20')


class D10Generator(bpy.types.Operator):
    """Generate a D10"""
    bl_idname = 'mesh.d10_add'
    bl_label = 'D10 Trapezohedron'
    bl_description = 'Generate an d10 trapezohedron dice'
    bl_options = {'REGISTER', 'UNDO'}

    size: FloatProperty(
        name='Face2face Length',
        description='Face-to-face size of the die',
        min=1,
        soft_min=1,
        max=100,
        soft_max=100,
        default=17,
        unit='LENGTH'
    )

    height: FloatProperty(
        name='Dice Height',
        description='Height of the die',
        min=0.45,
        soft_min=0.45,
        max=2,
        soft_max=2,
        default=2 / 3
    )

    add_numbers: BoolProperty(
        name='Generate Numbers',
        default=True
    )

    number_scale: FloatProperty(
        name='Number Scale',
        description='Size of the numbers on the die',
        min=0.1,
        soft_min=0.1,
        max=2,
        soft_max=2,
        default=1
    )

    number_depth: FloatProperty(
        name='Number Depth',
        description='Depth of the numbers on the die',
        min=0.1,
        soft_min=0.1,
        max=2,
        soft_max=2,
        default=0.75,
        unit='LENGTH'
    )

    font_path: StringProperty(
        name='Font',
        description='Number font',
        maxlen=1024,
        subtype='FILE_PATH'
    )

    number_v_offset: FloatProperty(
        name='Number V Offset',
        description='Vertical offset of the number positioning',
        min=0.0,
        soft_min=0.0,
        max=1,
        soft_max=1,
        default=1 / 3
    )

    number_indicator_type: EnumProperty(
        name='Orientation Indicator',
        items=(
            (NUMBER_IND_NONE, 'None', ','),
            (NUMBER_IND_BAR, 'Bar', ''),
            (NUMBER_IND_PERIOD, 'Period', ''),
        ),
        default=NUMBER_IND_PERIOD,
        description='Orientation indicator for numbers 6 and 9'
    )

    period_indicator_scale: FloatProperty(
        name='Period Scale',
        description='Scale of the period orientation indicator',
        min=0.1,
        soft_min=0.1,
        max=2,
        soft_max=2,
        default=1
    )

    period_indicator_space: FloatProperty(
        name='Period Space',
        description='Space between the period orientation indicator and the number',
        min=0,
        soft_min=0,
        max=3,
        soft_max=3,
        default=1
    )

    bar_indicator_height: FloatProperty(
        name='Bar Height',
        description='Height scale of the bar orientation indicator',
        min=0.1,
        soft_min=0.1,
        max=3,
        soft_max=3,
        default=1
    )

    bar_indicator_width: FloatProperty(
        name='Bar Width',
        description='Width scale of the bar orientation indicator',
        min=0.1,
        soft_min=0.1,
        max=2,
        soft_max=2,
        default=1
    )

    bar_indicator_space: FloatProperty(
        name='Bar Space',
        description='Space between the bar orientation indicator and the number',
        min=0,
        soft_min=0,
        max=3,
        soft_max=3,
        default=1
    )

    center_bar: BoolProperty(
        name='Center Align Bar',
        description='If true, the bar indicator is included in the vertical alignment of the number',
        default=True
    )

    one_offset: FloatProperty(
        name='Number 1 Offset',
        description='Offset the number 1 horizontally for an alternative centering',
        min=0,
        soft_min=0,
        max=1,
        soft_max=1,
        default=0
    )

    def execute(self, context):
        return execute_generator(self, context, D10Mesh, 'd10', height=self.height,
                                 number_v_offset=self.number_v_offset)


class D100Generator(bpy.types.Operator):
    """Generate a D100"""
    bl_idname = 'mesh.d100_add'
    bl_label = 'D100 Trapezohedron'
    bl_description = 'Generate an d100 trapezohedron dice'
    bl_options = {'REGISTER', 'UNDO'}

    number_indicator_type = NUMBER_IND_NONE
    one_offset = 0

    size: FloatProperty(
        name='Face2face Length',
        description='Face-to-face size of the die',
        min=1,
        soft_min=1,
        max=100,
        soft_max=100,
        default=17,
        unit='LENGTH'
    )

    height: FloatProperty(
        name='Dice Height',
        description='Height of the die',
        min=0.45,
        soft_min=0.45,
        max=2,
        soft_max=2,
        default=2 / 3
    )

    add_numbers: BoolProperty(
        name='Generate Numbers',
        default=True
    )

    number_scale: FloatProperty(
        name='Number Scale',
        description='Size of the numbers on the die',
        min=0.1,
        soft_min=0.1,
        max=2,
        soft_max=2,
        default=1
    )

    number_depth: FloatProperty(
        name='Number Depth',
        description='Depth of the numbers on the die',
        min=0.1,
        soft_min=0.1,
        max=2,
        soft_max=2,
        default=0.75,
        unit='LENGTH'
    )

    font_path: StringProperty(
        name='Font',
        description='Number font',
        maxlen=1024,
        subtype='FILE_PATH'
    )

    number_v_offset: FloatProperty(
        name='Number V Offset',
        description='Vertical offset of the number positioning',
        min=0.0,
        soft_min=0.0,
        max=1,
        soft_max=1,
        default=1 / 3
    )

    def execute(self, context):
        return execute_generator(self, context, D100Mesh, 'd100', height=self.height,
                                 number_v_offset=self.number_v_offset)


class MeshDiceAdd(Menu):
    """
    Dice menu under "Add Mesh"
    """

    bl_idname = 'VIEW3D_MT_mesh_dice_add'
    bl_label = 'Dice'

    def draw(self, context):
        layout = self.layout
        layout.operator_context = 'INVOKE_REGION_WIN'
        layout.operator('mesh.d4_add', text='D4 Tetrahedron')
        layout.operator('mesh.d4_crystal_add', text='D4 Crystal')
        layout.operator('mesh.d4_shard_add', text='D4 Shard')
        layout.operator('mesh.d6_add', text='D6 Cube')
        layout.operator('mesh.d8_add', text='D8 Octahedron')
        layout.operator('mesh.d10_add', text='D10 Trapezohedron')
        layout.operator('mesh.d100_add', text='D100 Trapezohedron')
        layout.operator('mesh.d12_add', text='D12 Dodecahedron')
        layout.operator('mesh.d20_add', text='D20 Icosahedron')


# Define "Extras" menu
def menu_func(self, context):
    layout = self.layout
    layout.operator_context = 'INVOKE_REGION_WIN'

    layout.separator()
    layout.menu('VIEW3D_MT_mesh_dice_add', text='Dice', icon='CUBE')


classes = [
    MeshDiceAdd,
    D4Generator,
    D4CrystalGenerator,
    D4ShardGenerator,
    D6Generator,
    D8Generator,
    D10Generator,
    D100Generator,
    D12Generator,
    D20Generator
]


def register():
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)

    # Add "Dice" menu to the "Add Mesh" menu
    bpy.types.VIEW3D_MT_mesh_add.append(menu_func)


def unregister():
    # Remove "Dice" menu from the "Add Mesh" menu.
    bpy.types.VIEW3D_MT_mesh_add.remove(menu_func)

    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(cls)
