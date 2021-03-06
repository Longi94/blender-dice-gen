# Blender Dice Gen

A Blender add-on that generates polyhedral dice.

![](https://github.com/Longi94/blender-dice-gen/raw/main/img/banner_ss.png)

## Installation and usage

Download the [python script](https://github.com/Longi94/blender-dice-gen/releases), then go to `Edit > Preferences... > Add-ons > Install...` and select the python file.

Afterwards you will be able to generate dice in `Add > Mesh > Dice`

The add-on will create 2 objects, a blank dice mesh and the numbers. The blank dice object will have a boolean modifier that has the "Realtime" flag turned off for performance reasons. Turn that on to see the result in the viewport.

Go to the [wiki](https://github.com/Longi94/blender-dice-gen/wiki/Properties) to see what each property does in detail.

## Supported dice

Create an issue if you'd like a new dice type to be added.

- D4 Tetrahedron
- D4 Crystal
- D4 Shard
- D6 Cube
- D8 Octahedron
- D10 Pentagonal Trepazohedron
- D100 Pentagonal Trepazohedron
- D12 Dodecahedron
- D20 Icosahedron

## Tips

- The unit scale in Blender is weird and with default settings the scale of the generated STL will be off by a factor of 1000 compared to the displayed scale in Blender. To have Blender display scales that matches the resulting STL set `Scene Properties > Units > Unit Scale` to 0.001. It also helps to set the length unit to millimeters.
- If the dice disappears when enabling the boolean modifier, try ticking the `Self` option under the solver options for the Exact solver. Or try switching to the Fast solver.

## Contributing

Please make pull requests against the `dev` branch to keep the `master` branch clean.
