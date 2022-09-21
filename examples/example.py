import numpy as np
import openmc

#INSTRUCTIONS: place these cross sections at the top of your OpenMC python script
#OR build your script from this baseline
#If you wish to reduce the number of angular moments, do so manually in the script
#OR alter the XS input file
####################################################################################################
###################################-Beginning of Cross Sections-####################################
####################################################################################################
#Group structure:
groups = openmc.mgxs.EnergyGroups([ 0.00000000E+00, 2.00000000E-02, 6.32455532E+02, 2.00000000E+07])
#Data for Material 1:
mat1_xsdat = openmc.XSdata('mat_1', groups)
mat1_xsdat.order = 4
mat1_xsdat.set_total([ 2.15783010E+01, 8.21277509E+01, 1.01310391E+03], temperature=294.)
mat1_xsdat.set_absorption([ 3.63281792E+00, 5.15434610E+01, 9.76506364E+02], temperature=294.)
mat1_xsdat.set_fission([ 3.20918209E+00, 3.65046490E+01, 8.34237159E+02], temperature=294.)
mat1_xsdat.set_nu_fission([ 8.23635331E+00, 8.89508781E+01, 2.03278569E+03], temperature=294.)
mat1_xsdat.set_chi([ 1.00000000E+00, 0.00000000E+00, 0.00000000E+00], temperature=294.)
scatter_matrix = np.array(\
    [[[ 1.79454827E+01, 0.00000000E+00, 0.00000000E+00],
      [ 3.58330074E-07, 3.04893077E+01, 0.00000000E+00],
      [ 0.00000000E+00, 9.49822669E-02, 3.65975415E+01]],
     [[ 5.46782575E+00, 0.00000000E+00, 0.00000000E+00],
      [-2.01191861E-07, 1.46335857E-01, 0.00000000E+00],
      [ 0.00000000E+00,-6.14166711E-02,-2.68010207E-01]],
     [[ 2.69534308E+00, 0.00000000E+00, 0.00000000E+00],
      [ 1.20232946E-09,-7.82543393E-02, 0.00000000E+00],
      [ 0.00000000E+00, 3.66341456E-02, 6.25168322E-02]],
     [[ 1.72286247E+00, 0.00000000E+00, 0.00000000E+00],
      [ 1.12560862E-07,-1.77310386E-01, 0.00000000E+00],
      [ 0.00000000E+00,-3.40446552E-02, 1.39513438E-01]],
     [[ 1.07042268E+00, 0.00000000E+00, 0.00000000E+00],
      [-9.98393014E-08, 4.73195156E-01, 0.00000000E+00],
      [ 0.00000000E+00, 2.62703141E-02,-7.80475470E-01]]])
scatter_matrix = np.transpose(scatter_matrix)
mat1_xsdat.set_scatter_matrix(scatter_matrix, temperature=294.)
#Data for Material 2:
mat2_xsdat = openmc.XSdata('mat_2', groups)
mat2_xsdat.order = 4
mat2_xsdat.set_total([ 6.43263071E+03, 1.23102693E+04, 1.72299270E+04], temperature=294.)
mat2_xsdat.set_absorption([ 1.06775052E+01, 2.14363767E+01, 1.69512795E+02], temperature=294.)
mat2_xsdat.set_fission([ 0.00000000E+00, 0.00000000E+00, 0.00000000E+00], temperature=294.)
mat2_xsdat.set_nu_fission([ 0.00000000E+00, 0.00000000E+00, 0.00000000E+00], temperature=294.)
mat2_xsdat.set_chi([ 1.00000000E+00, 0.00000000E+00, 0.00000000E+00], temperature=294.)
scatter_matrix = np.array(\
    [[[ 5.95079984E+03, 0.00000000E+00, 0.00000000E+00],
      [ 4.41369375E+02, 8.19167540E+03, 3.52281151E-01],
      [ 2.97839972E+01, 4.09715751E+03, 1.70600619E+04]],
     [[ 4.07530938E+03, 0.00000000E+00, 0.00000000E+00],
      [ 2.04970132E+02, 6.24630087E+03, 3.22514306E-01],
      [ 3.40152955E+00, 1.77280980E+03, 6.62091179E+03]],
     [[ 1.66774612E+03, 0.00000000E+00, 0.00000000E+00],
      [-4.43283976E+01, 3.44344502E+03, 2.70844211E-01],
      [-1.14211380E+01,-2.53136011E+02, 2.46523013E+03]],
     [[ 1.23899848E+02, 0.00000000E+00, 0.00000000E+00],
      [-1.18283958E+02, 1.11672364E+03, 2.09621510E-01],
      [-3.85242749E+00,-8.00041095E+02, 9.97129513E+02]],
     [[-2.35932234E+02, 0.00000000E+00, 0.00000000E+00],
      [-4.50908467E+01,-5.25990681E+01, 1.50806734E-01],
      [ 5.23484056E+00,-4.17300068E+02, 4.69602501E+02]]])
scatter_matrix = np.transpose(scatter_matrix)
mat2_xsdat.set_scatter_matrix(scatter_matrix, temperature=294.)
#Create the cross sections hdf5 file:
mg_cross_sections_file = openmc.MGXSLibrary(groups)
mg_cross_sections_file.add_xsdata(mat1_xsdat)
mg_cross_sections_file.add_xsdata(mat2_xsdat)
mg_cross_sections_file.export_to_hdf5('cross_sections.h5')
#Assign each cross section to a separate material:
materials = {}
for xs in ['mat_1','mat_2']:
    materials[xs] = openmc.Material(name=xs)
    materials[xs].set_density('macro', 1.)
    materials[xs].add_macroscopic(xs)
#Create the materials file for this specification:
materials_file = openmc.Materials(materials.values())
materials_file.cross_sections = 'cross_sections.h5'
materials_file.export_to_xml()
####################################################################################################
######################################-End of Cross Sections-#######################################
####################################################################################################
#The following creates the settings and sets the energy mode to multi-group:
settings_file = openmc.Settings()
settings_file.energy_mode = 'multi-group'
#The user should complete the OpenMC input below by specifying geometry, settings, etc.
