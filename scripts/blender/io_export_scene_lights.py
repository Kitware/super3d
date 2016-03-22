#!BPY
#ckwg +28
#Copyright 2014 by Kitware, Inc.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#  * Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
#  * Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
#  * Neither name of Kitware, Inc. nor the names of any contributors may be used
#    to endorse or promote products derived from this software without specific
#    prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ``AS IS''
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


bl_info = {
    'name': 'Export Scene Lights',
    'description': "Export scene light information (name, type, location, "
                   "rotations, scale)",
    'author': "Paul Tunison",
    'version': (0, 1),
    'blender': (2, 6, 9),
    'location': 'File > Export > Scene Lights',
    'category': 'Import-Export',
    }

import bpy
from bpy_extras.io_utils import ExportHelper
import mathutils


def export_scene_lights(bl_context, file_path):
    with open(file_path, 'w') as f:
        # Write out a header so the file makes some sense to a human
        f.write(
            "# Format:\n"
            "# +- lamp type\n"
            "# |   +- location       +- rotation quaternion  +- scale          +- look vector\n"
            "# |   |                 |                       |                 |\n"
            "# vvv vvvvvvvvvvvvvvvvv vvvvvvvvvvvvvvvvvvvvvvv vvvvvvvvvvvvvvvvv vvvvvvvvvvvvvvvvv\n"
            "# str float float float float float float float float float float float float float\n"
            "#\n"
            "# The look vector should be interpreted where the x-axis is positive towards\n"
            "# the top of the image, the y-axis to the left and the z-axis towards the\n"
            "# viewer.\n"
            "#\n"
            )

        lamps = []
        for obj in bl_context.scene.objects:
            if obj.type == 'LAMP':
                lamps.append(obj)
        print("Scene lamps:", lamps)
        for l in lamps:
            print('---', l.name, '---')

            obj_name = l.name                 # str
            print("\tobj name  :", obj_name)

            l_type = l.data.type            # str
            print("\tlamp type :", l_type)

            l_name = l.data.name       # str
            print("\tlamp name :", l_name)

            l_loc = l.location              # Vector
            print("\tlamp loc  :", l_loc)

            # force rotation mode to quaternion if not already, then switch
            # back, because quaternion value not computed until its the active
            # mode.
            orig_rot_mode = l.rotation_mode
            l.rotation_mode = 'QUATERNION'
            l_rot = l.rotation_quaternion   # Quaternion
            print("\tlamp quat :", l_rot)
            l.rotation_mode = orig_rot_mode

            l_scale = l.scale               # Vector
            print("\tlamp scale:", l_scale)

            l_look = mathutils.Vector()
            l_look.z = -1
            l_look.rotate(l_rot)

            f.write("%s %.14f %.14f %.14f %.14f %.14f %.14f %.14f %.14f %.14f %.14f %.14f %.14f %.14f\n"
                    % (l_type,
                       l_loc.x, l_loc.y, l_loc.z,
                       l_rot.w, l_rot.x, l_rot.y, l_rot.z,
                       l_scale.x, l_scale.y, l_scale.z,
                       l_look.x, l_look.y, l_look.z) )


class LightExporter (bpy.types.Operator):
    bl_idname = 'export_objects.scene_lights'
    bl_label = "Export Scene Lights"

    # "filepath" is apparently a magical property variable that, if defined, is
    # given the path to the file specified in the window... I cannot find how
    # this is actually set...
    filepath = bpy.props.StringProperty(subtype="FILE_PATH")

    def execute(self, context):
        export_scene_lights(context, self.filepath)
        return {'FINISHED'}

    def invoke(self, context, event):
        # Created a file selector window, which calls our execute function
        # upon loading.
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}


def menu_func(self, context):
    """ Call-back method for adding a menu option for this operation
    """
    self.layout.operator_context = 'INVOKE_DEFAULT'
    self.layout.operator(LightExporter.bl_idname, text="Scene Lights (.lights)")


###
# Registration
#
# (un)register_module(__name__) may also be used, but is more of a convenience
# approach when there are many classes being defined here so we don't have to
# make a (un)register_class call for each class defined. This will effectively
# cause blender to search this module for classes that sub-class registerable
# types and register them.
#
# As there in only one class being defined here, the lone (un)register_class
# calls are sufficient.
#

def register():
    bpy.utils.register_class(LightExporter)
    # adding menu option -- no explanation as to the use of the outside of the
    #                       example given in their docs.
    bpy.types.INFO_MT_file_export.append(menu_func)


def unregister():
    bpy.utils.unregister_class(LightExporter)
    # removing menu option -- no explanation as to the use of the outside of the
    #                         example given in their docs.
    bpy.types.INFO_MT_file_export.remove(menu_func)


# Only for user testing
if __name__ == '__main__':
  register()
