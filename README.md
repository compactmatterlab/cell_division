# cell_division


Model by Nadia Beydoun and Parag katira
Last Update 12/17/2020

Instructions

Download solver.m, dynamics.m, periphery.m

Load code into a folder

Select "Run". A movie and multiple plots from the data will be saved into the folder

The code determines the positions over time of a single cell rounding during mitosis using the ode 45 solver

Differential equations are in dynamics.m lines 324-345
         See differential equations to determine which forces affect the centrosomes, chromosomes, and periphery points
         Fc2p, Fcr2c and Frepc affect centrosomes
         Frepcr,Fcr2c, Frpa, and Frpe affect chromosomes
         Fatp, Fcph, Fcont, and F_pe affect periphery points. 


Each for-loop in dynamics.m represents a force on the elements 
         Each force is explained within the code. Force calculations can be manipulated based on future experimental research 
 

Any of the Boundary Conditions and Force Magnitudes along with constants in dynamics.m lines 104 and 86 can be manipulated 
         The function of these parameters are explained in the code and manipulation of these parameters affects their values in dynamics.m since they are all global variables.
         
         solver.m line 56 "pcell" is the main determinant for cellular volume 
         solver.m line 74 "fmagce" is the main determinant for centrosome separation


dynamics.m line 93 and 114 suppress z direction forces from coverplate and substrate on the total force calculation in line 117 in dynamics.  

Comment out dynamics.m line 112 when adding a cover plate to the cell. 

Uncomment solver.m line 518 if you want to caputure each frame of the movie

Uncomment solver.m line 521 if you want to plot the force due to pressure in the movie

Uncomment solver.m lines 555 to 623 if you want to plot the force vectors on the centrosome in the movie
         Comment out solver.m lines 585-592 and uncomment lines 594-597 to plot force vectors on periphery points in the movie




 
