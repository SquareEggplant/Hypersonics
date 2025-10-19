Documentation for 1-Dimensional Implicit Heat Analysis with Atmospheric Coupling
Overview	1
Program Structure	2
Function: altgen(inputfile)	2
Purpose	2
Inputs	2
Outputs	3
Process	3
Main Script: Heat Transfer Model	4
Purpose	4
Inputs	4
a. Material Properties (tungsten.txt)	4
b. Specific Heat Data (N2O.txt)	4
c. Atmospheric Data (thermal_team_data.csv)	4
Constants and Parameters	4
Numerical Implementation	5
a. Discretization	5
b. Boundary Conditions	5
c. Update Equations	5
d. Convective Heating Model	6
e. Radiative Heat Flux	6
Output	6
Variables	6
Visualization	6
Key Physical Assumptions	7






Overview
In response to the University Consortium for Applied Hypersonics’ (UCAH) call to design a hypersonic glide vehicle, the University of Maryland decided to assemble a team of undergraduates to create our own hypersonic vehicle. As per the solicitation sent by UCAH, the internal temperature requirement must not exceed 74°C (347 Kelvin). To better understand how commercial software like Ansys Thermal Desktop or COMSOL Multiphysics Software simulate heat transfer on 3-dimensional bodies, our team decided to create our own heat transfer solver. The purpose of this MATLAB script is to perform a 1-dimensional transient heat transfer simulation along a surface subjected to a high-speed flow. 
To conduct this analysis, we created an explicit Finite Difference Method (FDM). FDM solves a given partial differential equation (PDE) by estimating the gradients as a linear slope. We then apply a Neumann boundary condition at the leading edge (x = 0 meters) and a Dirichlet boundary condition at x = L meters, where L is the characteristic length of the body. At the Neumann boundary, we employ convection, radiation, and conduction with convection and conduction in the positive x-direction and radiation in the negative x-direction. At the Dirichlet boundary, we assume that the gradient of the temperature field is 0.
The convective heat transfer is modeled using Reynolds’ Analogy and aerothermal stagnation heating on a spherical surface. Reynolds’ Analogy estimates convective heat flux by computing the shear stress, 𝝉w, at the surface of the body. The advantage of using Reynolds’ Analogy is that the CFD software our group uses, Champs+, provides the shear stress at every point on the surface, so it is easy to implement into the design pipeline. However, the current model is meant for understanding the stagnation point heating, as that is where we will see our hottest temperatures. For this reason, we implement an equation that estimates the aerothermal heating at the stagnation point of a sphere. As for the radiative heat transfer, it propagates outward from the body as the heating due to radiation is much less than the radiation due to the temperature of the body. The conduction into the body is based on the thermal diffusivity which is calculated based on material properties, specified through a .txt file from Ansys Granta EduPack.






Program Structure
Trajectory Data
Generates atmospheric property fits (pressure, density, temperature, and Mach number) based on trajectory simulation results via a .csv file
Material Properties
Imports material values from Ansys Granta EduPack as a .txt file. This is required to calculate the conductive heat transfer, but is also a part of the stability condition.
Air Composition
Imports chemical data from NASA’s thermodynamic property data of common chemical compositions. Serves the purpose of estimating the specific heat capacity for constant pressure of the air.
Temperature Profile
Calculates the updated fit values at each time step and generates the temperature values at each point for every time step.
Assumptions
1D Heat Conduction:  Only solves for temperature gradients in the x-direction
Constant Material Properties: Material values do not vary with temperature currently
Fitting Values: Polynomial fits are not fully correct, estimate many of the values based on imperfect trajectory
Inaccurate Air Composition: Current data is based on N2O. There are many other chemical compounds that could be used but N2O is closest to the true composition of the air.
	




Function: altgen(inputfile)
Purpose
Computes the atmospheric pressure, density, temperature, and Mach number along the trajectory based on altitude data and fits them to high-order polynomial equations for smooth interpolation during simulation.
Inputs
inputfile (string): Path to a CSV file containing trajectory data with columns representing:


Time [s]


Altitude [m]


Mach number


(Unused placeholder)


(Unused placeholder)


Air density [kg/m³]


Air temperature [K]


Outputs
altitude — Matrix containing the original and computed atmospheric parameters.


rhofit — Polynomial coefficients for density vs. time.


tempfit — Polynomial coefficients for temperature vs. time.


localM — Polynomial coefficients for Mach number vs. time.


pfit — Polynomial coefficients for pressure vs. time.


Process
Reads and cleans the input table (readtable, rmmissing).


Converts data into numeric arrays.


Computes atmospheric pressure using the NASA Standard Atmosphere model:


Troposphere (0–11 km): temperature lapse rate.


Lower stratosphere (11–25 km): exponential decay.


Above 25 km: power-law relation.


Fits polynomial curves for each property using polyfit.



Main Script: Heat Transfer Model
Purpose
Simulates the temperature evolution of a 1D slab or nose tip surface under high-speed flight conditions, accounting for:
Conduction within the material.


Convective heating from aerodynamic flow.


Radiative cooling to the environment.



Inputs
a. Material Properties (tungsten.txt)
A text file in key-value format (e.g. Density_kg_m_3_: 19300) providing:
Density [kg/m³]


Thermal Conductivity [W/m·K]


Specific Heat Capacity [J/kg·K]


b. Specific Heat Data (N2O.txt)
A table of temperature vs. specific heat values, used to generate a polynomial fit for temperature-dependent Cp.
c. Atmospheric Data (thermal_team_data.csv)
Trajectory data passed to the altgen function.

Constants and Parameters
Variable
Description
Typical Value / Units
eps
Emissivity of Material
0.5
sig
Stefan–Boltzmann constant
5.67×10⁻⁸ W/m²·K⁴
L
Characteristic length
0.1 m
radius
Nose radius
0.003 m
Nx
Number of spatial nodes
1000
Fo
Fourier number (stability factor)
0.5
Pr
Prandtl number
0.71
r
Recovery factor
1.0
RgasAir
Gas constant for air
287.23 J/kg·K


Numerical Implementation
a. Discretization
The model divides the length L into Nx segments with spacing dx. Time is discretized using:
dt=dx2⋅Foαdt = \frac{dx^2 \cdot Fo}{\alpha}dt=αdx2⋅Fo​
where
α=kρcp\alpha = \frac{k}{\rho c_p}α=ρcp​k​
is the thermal diffusivity.
b. Boundary Conditions
Left boundary (x=0): Convective and radiative exchange with flow.
 Includes heat flux terms:
 qconv, qcond, qradq_{conv},\ q_{cond},\ q_{rad}qconv​, qcond​, qrad​
Right boundary (x=L): Insulated (Neumann condition),
 ∂T∂x=0\frac{\partial T}{\partial x} = 0∂x∂T​=0
c. Update Equations
Surface Node (x=0):
 Tnew(1)=T(1)+2 dtρcp(qconv−qrad)+2Fo(T(2)−T(1))T_{new}(1) = T(1) + \frac{2\,dt}{\rho c_p} (q_{conv} - q_{rad}) + 2Fo(T(2) - T(1))Tnew​(1)=T(1)+ρcp​2dt​(qconv​−qrad​)+2Fo(T(2)−T(1))
Internal Nodes:
 Tnew(i)=Fo[T(i+1)+T(i−1)]+(1−2Fo)T(i)T_{new}(i) = Fo[T(i+1) + T(i-1)] + (1-2Fo)T(i)Tnew​(i)=Fo[T(i+1)+T(i−1)]+(1−2Fo)T(i)
End Node (Insulated):
 Tnew(Nx)=T(Nx)+2Fo(T(Nx−1)−T(Nx))T_{new}(Nx) = T(Nx) + 2Fo(T(Nx-1) - T(Nx))Tnew​(Nx)=T(Nx)+2Fo(T(Nx−1)−T(Nx))
d. Convective Heating Model
Uses recovery temperature and stagnation properties derived from local Mach number and atmospheric fits:
Tstag=T∞(1+γ−12M2)T_{stag} = T_\infty \left(1 + \frac{\gamma - 1}{2} M^2\right)Tstag​=T∞​(1+2γ−1​M2)
Convective flux is then estimated as:
qconv=0.763 Pr−0.6(ρeμe)0.5(haw−hw)dudxq_{conv} = 0.763\,Pr^{-0.6}(\rho_e \mu_e)^{0.5}(h_{aw} - h_w)\sqrt{\frac{du}{dx}}qconv​=0.763Pr−0.6(ρe​μe​)0.5(haw​−hw​)dxdu​​
e. Radiative Heat Flux
qrad=εσT4q_{rad} = \varepsilon \sigma T^4qrad​=εσT4
Output
Variables
T_history: Matrix of temperatures [K] at each spatial node and time step.


T_new: Updated temperature profile per iteration.


Visualization
A 3D surface plot displays temperature evolution over time and position:
surf(x, t, T_history');
xlabel('Position [m]');
ylabel('Time [s]');
zlabel('Temperature [K]');
title('1D Transient Heat Transfer');
colorbar;
colormap hot;
shading interp;


