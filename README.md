# mo_adiabat

Compilation 
 - Please change the MPICXX in Makefile to your c++ compiler, and simply run make.

Running the program
- ./mosolid is the executable program. If you type in the potential temperature (K), the thermal profile will be calculated and shown in the order of T (K) / P (GPa). 
- n_pyrolite.txt (or a file set in main.cpp) sets the mantle composition. The composition is given in terms of the mol ratio (not wt%) of MgO-FeO-SiO2 ternary. The numbers do not need to add up to 100. It should work as long as the composition is reasonably close to pyrolite. I am not sure how it will behave when it is extremely enriched in FeO or SiO2.
- Please do not change the molecule_simple.txt. 

- The code can easily be expanded to other planets, but I haven’t made such changes… yet. 
