# MSci Project (in C++)


A group of C++ and associated header scripts, which encapsulate different approximations made for modelling the formation of giant planets in protoplanetary discs. (This is an updated version of [MSci_Project_Python](https://github.com/ocabrown/MSci_Project_Python), so the plots are very similar and I have removed some unnecessary scripts and classes).


## Data

The data the scripts use.


## Scripts

## C++ Scripts

### IsoT_Model.cpp (.hpp)
A class which allows a user to input initial conditions for the instantiation of a giant planet with an isothermal envelope. This simple model allowed us to test whether our code was working as expected as we can compare numerical outputs to their analytic equivalents.
### IsoD_Model.cpp (.hpp)
The same as IsoT_Model.py, but for an envelope with constant density.
### IsoO_Model.cpp (.hpp)
A class similar to above, but only uses a constant opacity assumption for the envelope. There are methods allowing users to run the equations governing the envelopes structure (Planet_Formation()), calculate the convergence stability of the code (L2norm()), and calculate the parameters which best describe the density profile (DensityFit()) - this would be used when determining the ideal initial conditions for the grid of mass values such that code stability is ensured.
### Freedman_Opacities_Function.cpp (.hpp)
A function to return the opacities calculated from the analytic fit from the Freedman et al. 2014 paper.
### FreO_Model.cpp (.hpp)
A class much like “IsoO_Model.cpp”, but using the Freedman opacity function to describe the opacity. There now is included a Brent’s method to find the temperature and opacity value for a given grid point simultaneously (as the opacity and temperature are now dependent upon one another).
### Combined_Opacities_Function.cpp (.hpp)
A script of functions which accumulates into the last function (Combined_Opacity()) which allows for a user to calculate an opacity from either the dust or gas depending on what temperature region the opacity is being calculated for (i.e. whether the dust has sublimed). The dust opacities are taken from Mordisini et al. 2012 paper b.
### ComO_Model.cpp (.hpp)
A class much like “FreO_Model.cpp”, but using the “Combined” opacity function to describe the opacity.
### Lum_Finder.cpp (.hpp)
A function which allows a user to find the luminosity value which causes the inner envelope radius to fall on a given giant planetary core. This luminosity can be used to determine how much mass the planet can accrete per unit time (thus evolve it).
### Evolution.cpp (.hpp)
A function which allows a user to evolve the envelope of a giant planet by increasing its mass according to its luminosity value and rerunning the respective opacity planet formation model’s method “Planet_Formation()”, and so on. This function can thus provide us with a mass accretion rate.
### Max_Acc_Rate_Calc.cpp (.hpp)
A function which calculates the maximum rate at which matter can be accreted by a forming giant planet, dependent upon how fast the surrounding disc can provide material. Different assumptions were made, informed from Tanigawa et al. 2012 and Lissauer et al. 2009, with the surface density profile of the disc informed by Crida et al. 2006.
### Mass_Grid_Solver.py
Functions that take the parameters output from DensityFit() methods in aforementioned classes and outputs a grid of mass values which would provide a logarithmically spaced grid of radius values for each mass value. (This would be repeatedly found until the mass grid doesn’t change - most stable grid using this method).
### Data_Output.cpp (.hpp)
A script of functions to output the data values from these prior scripts into .txt files on your Desktop. This data can then be imported into Python for plotting.
### Linspace.cpp (.hpp)
Functions to create vector of doubles that are spaced linearly and geometrically.
### Linear_Fit.cpp (.hpp)
Equivalent of Python’s scipy.optimize.curve_fit to calculate the parameters for DensityFit()
### Sign.cpp (.hpp)
Simple function to determine the sign of a double (1 or -1 for positive or negative)
### Interpolater.cpp (.hpp)
Interpolater function made in C++ (like Python’s interp1d).
### Convergence_Tester.cpp (.hpp)
A script with many functions to test the convergent properties of the different numerical integrations throughout this project.
### main.cpp
The script used to run all the previous classes and functions to output the data into a “Data_Outputs” file on desktop.

## Python Scripts (for plotting)

### Data_Input.py
Functions that extract and format data from external files into a more easily readable format for other Python scripts.
### Convergence_Plotter.py
A script that plots the output of Convergence_Tester.cpp.
### Max_Acc_Rate_Plotter.py
A script that plots the output of Max_Acc_Rate_Calc.cpp.
### L_Finder_Plotters.py
A script that plots the output of Lum_Finder.cpp.
### Class_Maker.py
A script that takes data inputs and forms classes of the data such that it can be used by Plotter.py.
### Tau_Line_Finder.py
A function to find where the optical depth is equal to 1 for a given envelope (used for plots).
### Plotter.py
A list of functions allowing a user to plot data outputted from the previous scripts (the same Plotter.py as is in the MSci_Project_Python).
### Run.py
A script to run all the inputting and plotting of the data to recreate the plots from MSci_Project_Python.
