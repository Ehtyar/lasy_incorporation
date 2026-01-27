# lasy_incorporation

## General

This repository contains three directories:
- bachelor_thesis contains the actual thesis and all the latex documents surrounding it as well as the images used and the references. It also contains the defense slides.
- python_modules contains all the python modules and jupyter notebooks developed for the bachelor thesis and afterwards fo the same project, as well as output folders for several of the modules.
- lib can be used as a python package containing all the functions developed for this work.

## Using the modules

To use the functions developed here, copy the lib directory or single files from it into your project directory or your python site-packages directory (and potentially rename it to fit your needs; lib is not a good name in most circumstances).

If you want to use showdata outside the package, you need to change the imports of .full_field and .ptime (which it needs) by removing the dot.

## Tests of the axiparabola laser

None of the tests done so far have successfully created an axiparabola flying focus laser, that actually moved as expected. However, there have been many smaller tests to figure out, why this is.

### In the Bachelors thesis

The bachelors thesis already contains a lot of tests for the modules in lib and the axiparabola laser in general. It can be found in bachlors_thesis/main.pdf. Here is an overview of the tests in it:
- A test of Lasy propagation using Gaussian pulses.
- A test of the full_field module for bringing Lasy laser pulses into PIConGPU simulations using Gaussian pulses.
- A test of the showdata module for generating the symlog plots seen in it and the functions show_w and plot_w for generating the w plots.
- Multiple tests of the RadialGroupDelay optical element found in the module of the same name, showing it does just delay parts of the pulse without focusing or defocusing.
- A test of the delay functions tau_D in the same module, showing they all line up, despite being described in different papers with different differential equations.
- A test showing the axiparabola does not work as expected.
- A test showing that it does not change anything to propagate the laser pulse to each point directly from the axiparabola instead of subsequently from om one point to the other.
- A test showing that it looks the same (or at least similarly not working) with PIConGPU simulations initialised with the laser already at the beginning of the focus region.
- A test showing, that the different equations for the shape of the axiparabola line upto one another and to the shape calculated by Lasy.
- A test showing, that this shapesolves the differential equation it is derived from.
- Multiple unsuccessful tests using the RGD and the axiparabola together.

### After the Bachelors thesis

The following test have been done afterwards. All the files and directories described here can be found in the python_modules directory.
- Multiple tests with different Lasy propagators, each with and without RGD:
  - AngularSpectrumPropagator:
    - Using the script runflfoc_angspec.py
    - Results in flfoc_angspec_out/ and rosi_flfoc_angspec_out/
    - Result: This propagator only works for very small numbers of points on the grid and is, therefore, unhelpful to this task
  - FresnelChirpZPropagator:
    - Using the script runflfoc_fresnel.py
    - Results in flfoc_fresnel_out/ and rosi_flfoc_fresnel_out/
    - In the images in flfoc_fresnel_out/ there are very weird lines visible next to the actually focused laser pulse.
    - Therefore, it is neccessary to use the grid_out option of the propagate function to only resample the field near the focus (see Lasy docu for information on how to do that).
    - The results of this, see rosi_flfoc_fresnel_out/, are very similar to the results of the standard propagator in Lasy (axiprop propagator). The main difference is the runtime per propagation step and because of this the number of steps.
- Because some papers mentioned using rectangular laser profiles (for their ease of calculation) I implemented the rectProfile module and tested with that:
  - The file flying_focus.ipynb (that was always the first test for new ideas and also contains other tests) shows, that the beam waist is larger with this new transverse profile (which is expected), but nothing else changes.
  - Testing more properly with the script flying_focus.py:
    - Results in flying_focus_img/
    - Although the beam waist is again larger than when testing with a transversal super-Gauss (as expected) the rest of the results looks very similar to what has been done earlier - no flying focus achieved.
- The beam waist in simulations with just the axiparabola looks similar to what one would expect from a parabolic mirror of similar focus length. Testing this:
  - The notebook read_the_file.ipynb shows the results of a PIConGPU simulation initialised with a laser pulse focused by just an axiparabola at the beginning of the focus region.
    - Near the end is a comparison between the beam waists in the simulation with the expectations of the axiparabola laser and a parabola mirror of the same focus length.
    - The measurement is more similar to the beam from the parabolic mirror, just with larger w0.
  - Testing this in Lasy using axiparabola_gauss.ipynb with a Gaussian pulse interacting with an axiparabola:
    - At first reproduced a very similar plot.
    - Then increased delta to make axiparabola and parabola more different - leads to the current state of the notebook.
    - The measured beam waist is now between the two descriptions (axiparabola theory and parabola with the same f0), however, it looks like it could be described by a parabola with a slightly longer f0 and a larger w0.
- Trying axiprop propagator with cartesian coordinates again but with more points, different starting beam profiles and again both with and without RGD:
  - Using the script runflfoc_axiprop_xyt.py.
  - Longitudinal Gauss, Transverse superGauss:
    - Results in rosi_sg_flfoc_axiprop_out/
    - (coming soon)
  - Longitudinal Gauss, Transverse rect:
    - Results in rosi_r1_flfoc_axiprop_out/
    - (coming soon)
  - Complete rect profile:
    - Results in rosi_r3_flfoc_axiprop_out/
    - (coming soon)
- Additional smaller tests to the ones mentioned in the Bachelors thesis can be found all over the python_modules directory in files and directories not directly mentioned in this section, as well as some additional images from tests in bachlors_thesis/Images/
  
### Open ideas

Some ideas for possible tests remain.
- Bringing the axiparabola near field directly into a PIConGPU simulation to see what will hapen - requires a large and pprobably very long simulation. Options:
  - Generate the pulse in Lasy, apply the axiparabola and save to openPMD compatible file and simulate in PIConGPU from there
  - Try to generate the pulse in PIConGPU analytically
- double checking the axiparabola derivation...
- trying to solve the propagation integrals analytically...
- idk
