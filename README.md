# Area-Preserving Parameterizations for Spherical Ellipses

[Ibón Guillén](https://webdiis.unizar.es/~iguillen/), [Carlos Ureña](https://lsi.ugr.es/curena/), Alan King, Marcos Fajardo, [Iliyan Georgiev](https://iliyan.com/), [Jorge López-Moreno](http://www.jorg3.com/), [Adrian Jarabo](https://webdiis.unizar.es/~ajarabo/).
In Computer Graphics Forum 36(4) 2017 (EGSR 2017)
[[Project page](http://giga.cps.unizar.es/~iguillen/projects/EGSR2017_Spherical_Ellipses/)]

## Overview
This repo contains an implementation of the disk solid angle techniques described in "Area-Preserving Parameterizations for Spherical Ellipses". The code is a fork of the Mitsuba 0.6.0 renderer (official repo: https://github.com/mitsuba-renderer/mitsuba).

The relevant files can be found at `src/shapes/disk`, which contains the different variants of the spherical ellipse sampling routines. It includes the precalculated solid angle tables for the tabulated variants.

## Compiling
This code has been tested on Ubuntu 16.04 and 18.04. Original Mitsuba code works on Windows and Mac too, so it should be possible to make it run on them.

  ### Linux (Ubuntu 16.04, 18.04)
   - `git clone https://github.com/iguillens/mitsuba-disk-sampling.git`
   - `sudo apt install git scons libboost-all-dev libpng-dev libjpeg-dev libopenexr-dev libxerces-c-dev libeigen3-dev libfftw3-dev libglewmx-dev freeglut3-dev`
   - `cd mitsuba-disk-sampling/`
   - `cp build/config-linux.py config.py`
   - `scons -j<number of cores>`

## Running
After a successful compilation, it should be enough to add `mitsuba` to the `$PATH` by
   - `source setpath.sh`

And run it
   - `mitsuba <path to scene>.xml`

## Example scenes
We include two of the scenes displayed in the publication (`spheres` and `bunny_media`). To test the different solid angle routines it is sufficient to change the `"sampling"` property of the disk area lights in the scenes with one of the following:
 - `"area"`  For regular area sampling.
 - `"gamito"`  For the rejection sampling method described in "Solid angle sampling of disk and cylinder lights" in Computer Graphics Forum 34(4) (2016).
 - `"polar"`  For the solid angle *polar* mapping variant.
 - `"polartab"`  For the tabulated solid angle *polar* mapping variant.
 - `"concentric"`  For the solid angle *concentric* mapping variant.
 - `"concentrictab"`  For the tabulated solid angle *concentric* mapping variant.
 - `"cylindrical"`  For the solid angle *parallel* mapping variant.
   
## Notes
 - By default the disk solid angle routines work in double precision, although it is possible to configure them to work in single precision in `src/shapes/disk.cpp`.
 - GUI is disabled by default to avoid compilation problems in modern systems.
  
If you have any questions and/or bug reports, feel free to contact Ibón Guillén (ibon@unizar.es).
