import math
import bpy
import os
from mathutils import Vector
from bpy.types import Menu
from bpy.props import FloatProperty, BoolProperty, StringProperty
from bpy_extras.object_utils import object_data_add

bl_info = {
    "name": "Dice Gen",
    "author": "Long Tran",
    "version": (0, 1, 0),
    "blender": (2, 80, 0),
    "location": "View3D > Add > Mesh",
    "description": "Generate dice meshes",
    "category": "Add Mesh",
}


def join(objects):
    bpy.context.view_layer.objects.active = objects[0]
    ctx = bpy.context.copy()
    # one of the objects to join
    ctx['active_object'] = objects[0]
    ctx['selected_editable_objects'] = objects
    bpy.ops.object.join(ctx)


def get_font(filepath):
    if filepath:
        bpy.ops.font.open(filepath=filepath)
        return next(filter(lambda x: x.filepath == filepath, bpy.data.fonts))
    else:
        return bpy.data.fonts['Bfont']


def create_number(context, number, font_path, size, number_depth, location, rotation):
    """
    Create a number mesh that will be used in a boolean modifier
    """
    # load the font
    font = get_font(font_path)

    # create the text curve
    font_curve = bpy.data.curves.new(type="FONT", name="number_{}".format(number))
    font_curve.body = str(number)
    font_curve.font = font
    font_curve.size = size

    # create object from curve
    curve_obj = bpy.data.objects.new("temp_number", font_curve)

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

    mesh_object.location.x = location[0] * size
    mesh_object.location.y = location[1] * size
    mesh_object.location.z = location[2] * size

    mesh_object.rotation_euler.x = rotation[0]
    mesh_object.rotation_euler.y = rotation[1]
    mesh_object.rotation_euler.z = rotation[2]

    # add solidify modifier
    bpy.ops.object.modifier_add(type='SOLIDIFY')
    context.object.modifiers[0].thickness = 2 * number_depth
    context.object.modifiers[0].offset = 0

    # hide object
    bpy.ops.object.hide_view_set(unselected=False)

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


class D6(bpy.types.Operator):
    """Generate a D6"""
    bl_idname = "mesh.d6_add"
    bl_label = "D6"
    bl_description = "Generate a cube dice"
    bl_options = {'REGISTER', 'UNDO'}

    size: FloatProperty(
        name="Size",
        description="Face-to-face size of the die",
        min=1,
        soft_min=1,
        max=100,
        soft_max=100,
        default=16,
        unit='LENGTH'
    )

    add_numbers: BoolProperty(
        name="Generate Numbers",
        default=True
    )

    number_depth: FloatProperty(
        name="Number Depth",
        description="Depth of the numbers on the die",
        min=0.1,
        soft_min=0.1,
        max=2,
        soft_max=2,
        default=0.75,
        unit='LENGTH'
    )

    font_path: StringProperty(
        name="Font",
        description="Number font",
        maxlen=1024,
        subtype='FILE_PATH',
    )

    def execute(self, context):
        # set font to emtpy if it's not a ttf file
        if self.font_path and os.path.splitext(self.font_path)[1] != '.ttf':
            self.font_path = ""
            return {'CANCELLED'}

        # create the cube mesh
        body_object = self.create_cube(context)

        # create number curves
        if self.add_numbers:
            self.create_numbers(context, body_object)

        return {'FINISHED'}

    def create_cube(self, context):
        # Calculate the necessary constants
        s = 0.5

        # create the vertices and faces
        v = [(-s, -s, -s), (s, -s, -s), (s, s, -s), (-s, s, -s), (-s, -s, s), (s, -s, s), (s, s, s), (-s, s, s)]
        faces = [[0, 3, 2, 1], [0, 1, 5, 4], [0, 4, 7, 3], [6, 5, 1, 2], [6, 2, 3, 7], [6, 7, 4, 5]]

        verts = [Vector(i) for i in v]

        # turn n-gons in quads and tri's
        faces = create_polys(faces)

        # resize
        verts = [i * self.size for i in verts]

        # generate object
        # Create new mesh
        mesh = bpy.data.meshes.new("Cube")

        # Make a mesh from a list of verts/edges/faces.
        mesh.from_pydata(verts, [], faces)

        # Update mesh geometry after adding stuff.
        mesh.update()

        return object_data_add(context, mesh, operator=None)

    def create_numbers(self, context, body_object):

        locations = [(0, 0, -0.5), (0, 0, 0.5)]
        rotation = [(0, math.pi, 0), (0, 0, 0)]

        numbers = []
        # create the number meshes
        for i in range(2):
            number_object = create_number(context, i + 1, self.font_path, self.size, self.number_depth, locations[i],
                                          rotation[i])
            numbers.append(number_object)

        # join the numbers into a single object
        join(numbers)
        numbers_object = context.view_layer.objects.active

        context.view_layer.objects.active = body_object

        # add boolean modifier
        bpy.ops.object.modifier_add(type='BOOLEAN')
        bpy.context.object.modifiers[0].object = bpy.data.objects[numbers_object.name]


class MeshDiceAdd(Menu):
    """
    Dice menu under "Add Mesh"
    """

    bl_idname = "VIEW3D_MT_mesh_dice_add"
    bl_label = "Dice"

    def draw(self, context):
        layout = self.layout
        layout.operator_context = 'INVOKE_REGION_WIN'
        layout.operator("mesh.d6_add", text="D6")


# Define "Extras" menu
def menu_func(self, context):
    layout = self.layout
    layout.operator_context = 'INVOKE_REGION_WIN'

    layout.separator()
    layout.menu("VIEW3D_MT_mesh_dice_add", text="Dice", icon="CUBE")


classes = [
    MeshDiceAdd,
    D6
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
