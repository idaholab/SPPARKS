# SPPARKS application RPV input  

  #seed	        ${SEED}  
  seed	        23810 

# define if 1NN or 2NN interaction 
  app_style       rpv 2 1 2 

  lattice         bcc 1.0
  region	    box block 0 40 0 40 0 4  
  # read_sites      100.in 

  create_box	box
  create_sites	box value i1 1
  set           i2 value 1 
  set           i2 value 2 if i1 = 1 fraction 0.2  # vacancy
  # set           i3 unique  
  # set           i2 value 3 if i1 = 1 fraction 0.0134 # Cu  
  # set         i2 value 4 if i1 = 1 fraction 0.01  # Ni
  # set         i2 value 5 if i1 = 1 fraction 0.015 # Mn 

  # sector          no
  sector          yes nstop 1 
  solve_style     tree 

# define 1NN and 2NN (optional) bond energy 
# ebond1 for 1NN; ebond2 for 2NN
# in the order of: 11 12 ... 1N; 21 22 ... 2N; ...; NN 

  ebond1          5 -0.611 -0.163 -0.480 -0.651 -0.446     0.126 -0.102 -0.213 -0.038    -0.414 -0.540 -0.365    -0.626 -0.366    -0.271    
  ebond2          5 -0.611 -0.163 -0.571 -0.596 -0.631    -0.014 -0.180 -0.193 -0.203    -0.611 -0.576 -0.621    -0.611 -0.496    -0.611 

#                 Fe        Cu        Ni        Mn 
  migbarrier      1 0.65    3 0.56    4 0.70    5 1.03
#  time_tracer     1.84  3 0.07 0.057  

# c11, c12, c44, ninteg, dcore 
# dtype (1 straight, 2 loop), burgers[3], xcore[3],line_vector,dradius,nsegment

#  moduli        243 145 116 21 0.36  
#  dislocation   1  1 0 0  0 0 0  1  5.0  40  

#  delta_v (eV/GPa)        V          Cu          Ni            Mn   (Si 0.00011) 
#  elastic_interaction   2 0.0  3 0.5   4  0.00075   5  0.00022  

#  sink_type, strength, shape, location[3], normal, radius, segment  
#  shape: 1 line; 2 polygon; 3 planar, 4 spherical 

#  sink            2  2.0  1  0 0 0  1  5.0  40  16.0

#  rsite rinput routput rbarrier rrate rtarget  

#  reaction        1  1 2  1.84  1.0e+8  10 

#  bfreq rpeak rdamp btypei btypej 

  ballistic   10 5.0 10.0 2 1 

# temperature in units of Kelvin  

  temperature	  873 

  diag_style      rpv stats yes list energy vac  
  stats           1000  
  dump            1 text 1000000 *.dump id i2 x y z d5 

  reset_time      0
  
  run             1000 
