# Blender Dice Gen

A Blender add-on that generates dice.

## Installation and usage

Download the python script, then go to `Edit > Preferences... > Add-ons > Install...` and select the python file.

Afterwards you will be able to generate dice in `Add > Mesh > Dice`

The add-on will create 2 objects, a blank dice mesh and the numbers. The blank dice object will have a boolean modifier that that has the "Realtime" flag turned off for performance reasons. Turn that on to see the result in the viewport.

## Tips

The unit scale in Blender is weird and with default settings the scale of the generated STL will be off by a factor of 1000 compared to the displayed scale in Blender. To have Blender display scales that matches the resulting STL set `Scene Properties > Units > Unit Scale` to 0.001. It also helps to set the length unit to millimeters.
