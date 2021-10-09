import math
import bpy
import os
from math import sqrt, acos
from mathutils import Vector, Matrix, Euler
from bpy.types import Menu
from bpy.props import FloatProperty, BoolProperty, StringProperty
from bpy_extras.object_utils import object_data_add
from add_mesh_extra_objects.add_mesh_solid import createSolid, createPolys

HALF_PI = math.pi / 2

# https://dmccooey.com/polyhedra
CONSTANTS = {
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
    }
}

bl_info = {
    'name': 'Dice Gen',
    'author': 'Long Tran',
    'version': (0, 1, 0),
    'blender': (2, 93, 0),
    'location': 'View3D > Add > Mesh',
    'description': 'Generate dice meshes',
    'category': 'Add Mesh',
}


class Mesh:

    def __init__(self, name):
        self.vertices = None
        self.faces = None
        self.name = name
        self.dice_mesh = None
        self.base_font_scale = 1

    def create(self, context):
        verts = [Vector(i) for i in self.vertices]

        # turn n-gons in quads and tri's
        faces = createPolys(self.faces)

        # generate object
        # Create new mesh
        mesh = bpy.data.meshes.new(self.name)

        # Make a mesh from a list of verts/edges/faces.
        mesh.from_pydata(verts, [], faces)

        # Update mesh geometry after adding stuff.
        mesh.update()

        self.dice_mesh = object_data_add(context, mesh, operator=None)
        return self.dice_mesh

    def get_numbers(self):
        raise []

    def get_number_locations(self):
        raise []

    def get_number_rotations(self):
        raise []

    def create_numbers(self, context, size, number_scale, number_depth, font_path):
        numbers = self.get_numbers()
        locations = self.get_number_locations()
        rotations = self.get_number_rotations()

        font_size = self.base_font_scale * size * number_scale

        numbers_object = create_numbers(context, numbers, locations, rotations, font_path, font_size, number_depth)

        apply_boolean_modifier(context, self.dice_mesh, numbers_object)


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
        return [str(i + 1) for i in range(6)]

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
        return [str(i + 1) for i in range(8)]

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
        return [str(i + 1) for i in range(12)]

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

        angles[5].y = -HALF_PI
        angles[5].rotate(Euler((-math.radians(72), 0, 0), 'XYZ'))
        angles[5].rotate(Euler((0, 0, (math.pi - CONSTANTS['octahedron']['dihedral_angle']) / 2), 'XYZ'))

        angles[6].y = -HALF_PI
        angles[6].rotate(Euler((math.radians(108), 0, 0), 'XYZ'))
        angles[6].rotate(Euler((0, 0, (math.pi - CONSTANTS['octahedron']['dihedral_angle']) / 2), 'XYZ'))

        angles[7].y = HALF_PI
        angles[7].rotate(Euler((math.radians(72), 0, 0), 'XYZ'))
        angles[7].rotate(Euler((0, 0, -(math.pi - CONSTANTS['octahedron']['dihedral_angle']) / 2), 'XYZ'))

        angles[8].z = math.radians(-36)
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
        return [str(i + 1) for i in range(20)]

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
    bpy.ops.object.join(ctx)


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


def apply_boolean_modifier(context, body_object, numbers_object):
    """
    Add a BOOLEAN modifier to body_object that targets
    :param context:
    :param body_object:
    :param numbers_object
    :return:
    """
    context.view_layer.objects.active = body_object
    bpy.ops.object.modifier_add(type='BOOLEAN')
    bpy.context.object.modifiers[0].object = bpy.data.objects[numbers_object.name]
    bpy.context.object.modifiers[0].show_viewport = False


def create_numbers(context, numbers, locations, rotations, font_path, font_size, number_depth):
    number_objs = []
    # create the number meshes
    for i in range(len(locations)):
        number_object = create_number(context, numbers[i], font_path, font_size, number_depth, locations[i],
                                      rotations[i])
        number_objs.append(number_object)

    # join the numbers into a single object
    join(numbers)
    apply_transform(context.view_layer.objects.active, use_rotation=True)
    return context.view_layer.objects.active


def create_number(context, number, font_path, font_size, number_depth, location, rotation):
    """
    Create a number mesh that will be used in a boolean modifier
    """
    # load the font
    font = get_font(font_path)

    # create the text curve
    font_curve = bpy.data.curves.new(type='FONT', name=f'number_{number}')
    font_curve.body = str(number)
    font_curve.font = font
    font_curve.size = font_size

    # create object from curve
    curve_obj = bpy.data.objects.new('temp_number', font_curve)

    # convert curve to mesh
    mesh = curve_obj.to_mesh().copy()
    curve_obj.to_mesh_clear()
    bpy.data.objects.remove(curve_obj)
    bpy.data.curves.remove(font_curve)

    # add number
    mesh_object = object_data_add(context, mesh, operator=None)

    # set origin to bounding box center
    mesh_object.select_set(True)
    bpy.ops.object.origin_set(type='ORIGIN_GEOMETRY', center='BOUNDS')

    mesh_object.location.x = location[0]
    mesh_object.location.y = location[1]
    mesh_object.location.z = location[2]

    mesh_object.rotation_euler.x = rotation[0]
    mesh_object.rotation_euler.y = rotation[1]
    mesh_object.rotation_euler.z = rotation[2]

    # add solidify modifier
    bpy.ops.object.modifier_add(type='SOLIDIFY')
    context.object.modifiers[0].thickness = 2 * number_depth
    context.object.modifiers[0].offset = 0

    return mesh_object


class D6Generator(bpy.types.Operator):
    """Generate a D6"""
    bl_idname = 'mesh.d6_add'
    bl_label = 'D6'
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
        default=1,
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
        subtype='FILE_PATH',
    )

    def execute(self, context):
        # set font to emtpy if it's not a ttf file
        self.font_path = validate_font_path(self.font_path)

        # create the cube mesh
        die = Cube('d6', self.size)
        die.create(context)

        # create number curves
        if self.add_numbers:
            die.create_numbers(context, self.size, self.number_scale, self.number_depth, self.font_path)

        return {'FINISHED'}


class D8Generator(bpy.types.Operator):
    """Generate a D8"""
    bl_idname = 'mesh.d8_add'
    bl_label = 'D8'
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
        default=1,
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
        subtype='FILE_PATH',
    )

    def execute(self, context):
        # set font to emtpy if it's not a ttf file
        self.font_path = validate_font_path(self.font_path)

        # create the cube mesh
        die = Octahedron('d8', self.size)
        die.create(context)

        # create number curves
        if self.add_numbers:
            die.create_numbers(context, self.size, self.number_scale, self.number_depth, self.font_path)

        return {'FINISHED'}


class D12Generator(bpy.types.Operator):
    """Generate a D12"""
    bl_idname = 'mesh.d12_add'
    bl_label = 'D12'
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
        default=1,
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
        subtype='FILE_PATH',
    )

    def execute(self, context):
        # set font to emtpy if it's not a ttf file
        self.font_path = validate_font_path(self.font_path)

        # create the cube mesh
        die = Dodecahedron('d12', self.size)
        die.create(context)

        # create number curves
        if self.add_numbers:
            die.create_numbers(context, self.size, self.number_scale, self.number_depth, self.font_path)

        return {'FINISHED'}


class D20Generator(bpy.types.Operator):
    """Generate a D20"""
    bl_idname = 'mesh.d20_add'
    bl_label = 'D20'
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
        default=1,
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
        subtype='FILE_PATH',
    )

    def execute(self, context):
        # set font to emtpy if it's not a ttf file
        self.font_path = validate_font_path(self.font_path)

        # create the cube mesh
        die = Icosahedron('d20', self.size)
        die.create(context)

        # create number curves
        if self.add_numbers:
            die.create_numbers(context, self.size, self.number_scale, self.number_depth, self.font_path)

        return {'FINISHED'}


class MeshDiceAdd(Menu):
    """
    Dice menu under "Add Mesh"
    """

    bl_idname = 'VIEW3D_MT_mesh_dice_add'
    bl_label = 'Dice'

    def draw(self, context):
        layout = self.layout
        layout.operator_context = 'INVOKE_REGION_WIN'
        layout.operator('mesh.d6_add', text='D6 Cube')
        layout.operator('mesh.d8_add', text='D8 Octahedron')
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
    D6Generator,
    D8Generator,
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
