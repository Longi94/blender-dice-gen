import math
import bpy
import os
from mathutils import Vector, Matrix
from bpy.types import Menu
from bpy.props import FloatProperty, BoolProperty, StringProperty
from bpy_extras.object_utils import object_data_add

HALF_PI = math.pi / 2

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
        faces = create_polys(self.faces)

        # generate object
        # Create new mesh
        mesh = bpy.data.meshes.new(self.name)

        # Make a mesh from a list of verts/edges/faces.
        mesh.from_pydata(verts, [], faces)

        # Update mesh geometry after adding stuff.
        mesh.update()

        self.dice_mesh = object_data_add(context, mesh, operator=None)
        return self.dice_mesh

    def get_number_locations(self):
        raise NotImplementedError()

    def get_number_rotations(self):
        raise NotImplementedError()

    def create_numbers(self, context, size, number_scale, number_depth, font_path):
        locations = self.get_number_locations()
        rotations = self.get_number_rotations()

        font_size = self.base_font_scale * size * number_scale

        numbers_object = create_numbers(context, locations, rotations, font_path, font_size, number_depth)

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
        self.v_coord_const = (size * math.sqrt(3)) / 2
        s = self.v_coord_const

        # create the vertices and faces
        self.vertices = [(s, 0, 0), (-s, 0, 0), (0, s, 0), (0, -s, 0), (0, 0, s), (0, 0, -s)]
        self.faces = [[4, 0, 2], [4, 2, 1], [4, 1, 3], [4, 3, 0], [5, 2, 0], [5, 1, 2], [5, 3, 1], [5, 0, 3]]

        self.base_font_scale = 0.7

    def get_number_locations(self):
        c = self.v_coord_const / 3
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


def create_numbers(context, locations, rotations, font_path, font_size, number_depth):
    numbers = []
    # create the number meshes
    for i in range(len(locations)):
        number_object = create_number(context, i + 1, font_path, font_size, number_depth, locations[i],
                                      rotations[i])
        numbers.append(number_object)

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


def create_polys(poly):
    """
    this function creates a chain of quads and, when necessary, a remaining tri
    for each polygon created in this script. be aware though, that this function
    assumes each polygon is convex.

    :param poly: list of faces, or a single face, like those needed for mesh.from_pydata.
    :return: the tessellated faces.
    """
    # check for faces
    if len(poly) == 0:
        return []
    # one or more faces
    if type(poly[0]) == type(1):
        poly = [poly]  # if only one,  make it a list of one face
    faces = []
    for i in poly:
        L = len(i)
        # let all faces of 3 or 4 verts be
        if L < 5:
            faces.append(i)
        # split all polygons in half and bridge the two halves
        else:
            f = [[i[x], i[x + 1], i[L - 2 - x], i[L - 1 - x]] for x in range(L // 2 - 1)]
            faces.extend(f)
            if L & 1 == 1:
                faces.append([i[L // 2 - 1 + x] for x in [0, 1, 2]])
    return faces


class D6Generator(bpy.types.Operator):
    """Generate a D6"""
    bl_idname = 'mesh.d6_add'
    bl_label = 'D6'
    bl_description = 'Generate a cube dice'
    bl_options = {'REGISTER', 'UNDO'}

    base_font_scale = 1

    size: FloatProperty(
        name='Size',
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
        cube = Cube('d6', self.size)
        cube.create(context)

        # create number curves
        if self.add_numbers:
            cube.create_numbers(context, self.size, self.number_scale, self.number_depth, self.font_path)

        return {'FINISHED'}


class D8Generator(bpy.types.Operator):
    """Generate a D8"""
    bl_idname = 'mesh.d8_add'
    bl_label = 'D8'
    bl_description = 'Generate a octahedron dice'
    bl_options = {'REGISTER', 'UNDO'}

    base_font_scale = 1

    size: FloatProperty(
        name='Size',
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
        octahedron = Octahedron('d8', self.size)
        octahedron.create(context)

        # create number curves
        if self.add_numbers:
            octahedron.create_numbers(context, self.size, self.number_scale, self.number_depth, self.font_path)

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
        layout.operator('mesh.d6_add', text='D6')
        layout.operator('mesh.d8_add', text='D8')


# Define "Extras" menu
def menu_func(self, context):
    layout = self.layout
    layout.operator_context = 'INVOKE_REGION_WIN'

    layout.separator()
    layout.menu('VIEW3D_MT_mesh_dice_add', text='Dice', icon='CUBE')


classes = [
    MeshDiceAdd,
    D6Generator,
    D8Generator
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
