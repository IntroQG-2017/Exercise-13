# Exercise 13: Quantitative thermochronology, part I
This exercise is part 1 of the exercises on thermochronology.
In this exercise you will run a time-dependent 1D heat transfer model, plot predicted geotherms, and predict thermochronometer ages for several thermochronometers.

## Overview
This two-part set of exercises is designed to give you a better understanding of thermochronology and the temperature field in the Earth's crust.
Thermochronology combines many aspects of what we have studied in this course, including the advection and diffusion equations, erosional processes, and some basic geostatistics.
This first part of the exercises on thermochronology is intended to provide an understanding of heat advection and diffusion in the Earth's crust, the time-dependence of thermal processes, and how the thermal history recorded in our "digital rock" can be used to predict thermochronometer ages.

## Time-dependent temperatures in the Earth
In this exercise we will use an analytical solution to the 1-D time-dependent thermal advection-diffusion equation to simulate erosion of rock at the Earth's surface, the upward transport of the underlying rock toward the surface, and the changes in a 1-D geotherm with time.
In the [starter code](age_predict_1D.py) for this week's exercise, we provide a working solution to the basic equation for calculating temperature *T* as a function of depth *z* and time *t* (*you don't need to convert the equation to Python, don't worry*).
The equation was originally published by Carslaw and Jaeger (1959)

![Equation 1](Images/Equation1.png)<br/> *Equation 1. 1D time-dependent heat advection-diffusion equation.*

where *G* is the initial temperature gradient (increase in temperature with depth), *v*<sub>*z*</sub> is the vertical advection velocity (positive upward), *κ* is the thermal diffusivity and `erfc()` is the complementary error function, defined as

![Equation 2](Images/Equation2.png)<br/> *Equation 2. The complementary error function.*

With this equation, the temperature increases linearly with depth at the start of the calculation (*T*(*z*,*t*=0) = *Gz*) and the geotherm will evolve over time as a function of the advection velocity and thermal diffusivity.

If you are curious about how the the advection and diffusion equations are combined in Equation 1, you can find a solution to the steady-state version of the 1D heat transfer equation in the [notes on solving the advection-diffusion equation](https://introqg.github.io/2017/lessons/L11/solving-advection.html) from [Lesson 11](https://introqg.github.io/2017/lessons/L11/overview.html).

## Problem 1 - The time dependence of rock advection
We will begin by plotting geotherms to get a sense of how our temperature equation works.

1. If you download [a copy of the Python script `age_predict_1D.py`](age_predict_1D.py) you should be able to run it without making any changes to produce a plot like that shown below.

    ![Time-dependent heat transfer with advection](Images/1D_transient_plot1.png)<br/>
    *Figure 1. 1D transient thermal solution including advection.*<br/><br/>
To start, please add axis labels and a title to this plot. **What is the advection velocity for this temperature calculation?** Add a text label listing the advection velocity on your the plot using the `plt.text()` function, then save a copy of the plot at the end of this document.
2. **How does the thermal solution change when you change the advection velocity?** Increase the advection velocity to 1.0 mm/a and save the resulting plot at the end of this document. Do the same thing for an advection velocity of 0.1 mm/a. **How does changing the advection velocity change temperatures in the shallow crust (<10 km depth)?**
3. Reset the advection velocity to 0.5 mm/a. Now increase the total simulation time to 100 Ma. Once again, save a copy of this plot at the end of this document. **What happens to temperatures in the model as the simulation time increases?** You might want to run some additional calculations (you don't need to save the plots) with even longer simulation times (1000 Ma, 5000 Ma, etc.). **What do you observe for the temperatures in the model? Are there any potential problems with these temperatures? Does the model approach a thermal steady state?**
4. Reset the simulation time to 50 Ma. Now we'll explore the effect of the initial temperature gradient. Increase the initial temperature gradient to 20°C/km, run the model and save the plot. Do the same for a temperature gradient of 5°C/km. **How does the initial temperature gradient affect the temperatures in the model? Is there any clear relationship between the initial temperature gradient and the maximum temperature in the model at *t* = 0 Ma?**
5. Lastly, reset the initial temperature gradient to 10°C/km. Increase the number of temperature calculations to plot from 1 to 5 and run the thermal model. Save this plot and insert it at the end of this document. **In your opinion, is it helpful to see the temperature calculations at different times?**

## Predicting thermochronometer ages
Thermochronometers record the time since a rock or mineral was at a certain temperature (the closure temperature) in the Earth.
Above, we have calculated temperature solutions assuming 1-D vertical advection of the crust.
Here, we will track a parcel of rock through the temperature field as it is advected toward the surface and use the recorded temperature history to predict thermochronometer ages.
The key to understanding what is done here is to understand that we will be simulating the position of a parcel of rock at depth in the earth, and at each time step the position of the rock parcel will be moved upward according to the length of the time step multiplied by the advection velocity.
In mathematical terms, this relationship is

![Equation 3](Images/Equation3.png)<br/> *Equation 3. Equation for calculating rock particle depth as a function of time.*

Thus, we will track a parcel of rock from some depth in the model at the start of the temperature calculation to its final location at the surface when the simulation is complete at 0 Ma.
At each depth, the temperature of the parcel of rock will be stored, which will allow the cooling in the temperature history of the rock parcel to be used to predict different thermochronometer ages.
We will consider the apatite (U-Th)/He, zircon (U-Th)/He, and muscovite <sup>40</sup>Ar/<sup>39</sup>Ar thermochronometers that were presented briefly in lecture.

Thermochronometer closure temperatures will be predicted using Dodson's method, which was also discussed briefly in lecture.
According to Dodson's method, the closure temperature *T*<sub>c</sub> of a thermochronometer is

![Equation 4](Images/Equation4.png)<br/> *Equation 4. The effective closure temperature according to Dodson's method.*

where *E*<sub>a</sub> is the activation energy, *R* is the universal gas constant, *A* is a geometric factor (*A* = 25 for a sphere, *A* = 8.7 for a planar sheet), *τ* is time for the diffusivity to decrease by a factor of 1/e, *D*<sub>0</sub> is the diffusivity at infinite temperature and *a* is the diffusion domain (we'll assume this is the size of the mineral). The value of *τ* can be calculated as a function of the cooling rate *dT*/*dt*

![Equation 5](Images/Equation5.png)<br/> *Equation 5. The characteristic time for a change in diffusivity.*

By simulating cooling of the minerals by iterating over the values of the recorded temperature history from depth to the surface, thermochronometer ages can be predicted for various systems.

## Problem 2 - Cooling ages and their relationship exhumation rates
One of the main interests for geoscientists using thermochronology is to determine the average exhumation rate of rocks in a study area based on the ages of thermochronometer data at the surface.

1. At the top of the Python script [`age_predict_1D.py`](age_predict_1D.py), there are three flags (`True`/`False` variables) that allow you to enable the calculation of different thermochronometers.
Set the value for `calc_AHe` to `True` to enable prediction of apatite (U-Th)/He ages and run the model with the default parameters.
**What is the predicted apatite (U-Th)/He age?**
**What does this age mean?**
**What is the closure temperature for the model in this case?**
If you look carefully through the code, you can find where the closure temperature is calculated.
Add text to your plot to display the predicted apatite (U-Th)/He age and predicted closure temperature using the `plt.text()` function.
I suggest that you use the plotting text function at the bottom of the script where the predicted apatite (U-Th)/He age is written to the screen (if requested).
Save a copy of the plot and insert it at the end of this document.
2. Similar to above, set the flags for `calc_ZHe` and `calc_MAr` to `True` and add the corresponding `plt.text()` functions to display the predicted zircon (U-Th)/He and muscovite <sup>40</sup>Ar/<sup>39</sup>Ar ages and their predicted closure temperatures.
As above, I suggest that you use the plotting text function at the bottom of the script where the predicted ages are written to the screen (if requested).
Save a copy of the plot and insert it at the end of this document.
**Do all of the predicted thermochronometer ages make sense?**
3. To understand more about how the thermochronometer ages are predicted, it can he helpful to look at the temperature-depth history of the parcel of rock as it travels from depth to the surface.
You can see this history visually by setting the flag `plot_Tzhist` to `True` at the top of the Python script.
This will add a set of symbols to the plot that display how the temperature and depth of the rock parcel changes with time.
**Considering the predicted thermochronometer ages and closure temperatures, do you see any problems with the predicted ages now?**
Save a copy of this plot and insert it at the end of this document.
4. Finally, increase the advection velocity in the thermal model to 1.0 mm/a and produce a similar plot to above.
**Do the predicted thermochronometer ages make more sense now?**
**What was the problem in questions 2 and 3?**
Save a copy of this plot and insert it at the end of this document.

## What to submit
**For this exercise, your modifications to the end of this document should include**

1. The 7 plots requested for problem 1 and the 4 plots requested for problem 2.
2. Figure captions for each plot describing the plot as if it were in a scientific journal article.
3. Answers to all of the questions in bold
4. Copies of your modified Python scripts for Problems 1 and 2

**NOTE**: You may want to reference these plots in your final report on these exercises, so be sure to keep the copies of the plot files.

## References
Carslaw, H. S., & Jaeger, J. C. (1959). Conduction of heat in solids. Oxford: Clarendon Press.

Dodson, M. H. (1973). Closure temperature in cooling geochronological and petrological systems. Contributions to Mineralogy and Petrology, 40(3), 259–274.

# Answers
## Problem 1
This is some text. You can use *italics* or **bold** text easily. You may want to read a bit more about [formatting text in Github-flavored Markdown](https://help.github.com/articles/basic-writing-and-formatting-syntax/). You can see an example of how to display an image with a caption below.

![Text shown if image does not load](Images/sine.png)<br/>
*Figure 2: Sine wave calculated from 0 to 2π*
