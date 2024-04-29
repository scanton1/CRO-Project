# CRO-Project
Work from Simon Anton for the CRO desalination project. 

## Goal
Give an overview of how the code and simulations work. The document will first discuss each UDF in the "cro.c" file and then how to run the simulations.

![image](https://github.com/scanton1/CRO-Project/assets/168468147/bf051a82-b992-40f2-82a3-327b435049a6)

# [User Defined Functions (UDFs)](https://www.afs.enea.it/project/neptunius/docs/fluent/html/udf/main_pre.htm)
All the UDFs necessary for running the simulations are found in "cro.c". I'll give a brief overview of the purpose and methodology of each function.

## DEFINE_ON_DEMAND(species_list)
This is a simple utility function that prints out each species and the corresponding integer index (i) in a mixture defined in Fluent. This function can be executed at any point manually by the user using the Fluent GUI. This UDF is useful to make sure you are indexing the right species when using certain UDF macros. For example, if you wanted species mass fraction cell data, you would use the macro C_YI(c,t,i). The index i needs to correspond to the right species in your mixture template. Thus, you would use this DEFINE_ON_DEMAND UDF to check the index associated with the desired species.

## DEFINE_ON_DEMAND(face_conc)
This is another utility tool that is simply used for printing data in a cell or face thread. The function loops through each element in that thread and prints whatever you desire. I used this UDF and alter it frequently to do sanity checks. 

## DEFINE_INIT(init_conc, d)
This UDF is used to initialize the species mass fraction in the computational domain. In the channel simulations, there were just two major regions that we needed to initialize: the feed channel and the membrane / porous zone. This UDF loops through all the cell threads in the domain and only initializes a non-zero value for the zones contained within the feed. 

## DEFINE_PROFILE(outlet_flux, t, yi)
This UDF is assigned to the feed outlet and approximates zero diffusive flux for the salt species ($\frac{\partial{C}}{\partial{x}}=0$). The UDF loops over all faces on the outlet and finds the adjacent cell to that face. Zero diffusive flux is then approximated by setting the face mass fraction to be equal to the cell mass fraction of the salt species.  

## DEFINE_PROPERTY(density, c, t)
This UDF is used when defining the mixture template to provide an expression for density as a function of salt concentration. This equation comes from a paper using sodium sulfate as the salt species. This is not the case for future simulations, so this needs to be updated for our ocean water solution. 

## DEFINE_DIFFUSIVITY(diffusivity, c, t, i)
This UDF is used when defining the mixture template to provide an expression for diffusivity as a function of salt concentration. This equation comes from a paper using sodium sulfate as the salt species. This is not the case for future simulations, so this needs to be updated for our ocean water solution.

## DEFINE_PROPERTY(viscosity, c, t)
This UDF is used when defining the mixture template to provide an expression for viscosity as a function of salt concentration. This equation comes from a paper using sodium sulfate as the salt species. This is not the case for future simulations, so this needs to be updated for our ocean water solution.

## DEFINE_PROFILE(porous_zone, t, i)
This UDF is used to make the permeability of the membrane a function of viscosity (which is in itself a function of concentration). The expression `1/(K*C_MU_EFF(c,t)*dz)` is derived from the paper referenced. This will likely change when simulating the membrane we bought. 

## DEFINE_PROFILE(membrane_concentration, t, y_i)
This UDF is one of the two fundamental functions to modeling RO desalination. This function takes the equation $JC_{m}+D\frac{\partial{C}}{\partial{x}}=JC_{p}$, discretizes the gradient term and solves for the concentration at the membrane $C_{m}$: $C_{m}=\frac{D(4C_{A}-C_{B})}{3D-2J\Delta{z}}$. Unfortunately, there is not a way to apply a DEFINE_PROFILE UDF to the faces of the membrane on an interior face zone. So, we elect to create an explicit cell zone directly over the membrane (membrane adjacent zone) to apply our profile UDF to. This may not be the best way, but we haven't not been able to find a different way to make it work.

For the $C_{m}$ equation above, $C_{A}$ and $C_{B}$ correspond to cells in two explicit cell zones titled cell zone A and cell zone B, respectively. These explicit cell zones have to be created while constructing your computational domain. This will be covered later in the simulations section. The whole purpose of this UDF is to assign the salt mass fraction of salt for a cell next to the membrane to the expression above for $C_{m}$. To do this you must find the values for $C_{A}$ and $C_{B}$. There is no easy way of doing this with built-in UDF macros, so this code finds those values then assigns $C_{m}$.

The function first defines the cell threads that correspond to cell zone A and cell zone B. The Fluent GUI allows you to find the zone ids for those zones easily. From there you can utilize a macro to assign a pointer to the cell thread. The function then loops over all cells in the cell zone adjacent to the membrane. For each cell in the membrane adjacent zone, the function then creates two separate loops over cell zone A and cell zone B to find the nearest cell in each respective zone from the current cell in the membrane adjacent loop. The closest cell from each region is what is referenced in the equation for $C_{m}$. Once $C_{A}$ and $C_{B}$ are found, we assign $C_{m}$. 

## DEFINE_PROFILE(permeate_pressure, t, p)
This UDF is the second fundamental function for modeling RO desalination. This function calculates the osmotic pressure and assigns it to the pressure outlet at the backside of the membrane. Although unclear from this, this function is what enforces the permeate flux through the membrane. This is due to how Fluent handles a [porous media](https://www.afs.enea.it/project/neptunius/docs/fluent/html/ug/node233.htm). For laminar flows through porous media, the pressure drop is simply Darcy's Law. Therefore, to get the desired flow rate through the membrane, all we need to do is control the pressure differential. The image below hopefully helps you see the equivalence between the flux equation $J=K(\Delta{p}-\Delta{\Pi})$ and the porous media treatment of pressure drop. The membrane constant $K$ can be expressed as $K=\frac{\alpha}{\mu \Delta{n}}$ where $\alpha$ is the permeability, $\mu$ is the viscosity, and $\Delta n$ is the thickness of the membrane.  

![image](https://github.com/scanton1/CRO-Project/assets/168468147/24f7c0cf-36b0-4bcf-a552-f6423475253d)

All the function does is loop through the faces on the outlet and find the nearest cell in the membrane adjacent cell zone for each of these faces, then it calculates the osmotic pressure difference assuming the permeate concentration is zero. This osmotic pressure difference is then assigned to the outlet face. 

## DEFINE_EXECUTE_AT_END(vol)
This UDF is executed at the end of each time step to calculate the volume of permeate produced over that time step. The function simply loops over the faces on the membrane outlet and calculates the volume "produced" by each cell using the velocity, face area and time step. This volume is then written to a file that can be summed in post processing to find the total volume produced. 

### Note
The functions above work but need to be altered depending on the zone ids for the different regions in your domain. This code should work in 3D as is, but it needs to be updated for our specific RO membrane and our salt water solution. 

# Simulation Setup
To recreate my simulations I will go over the three steps in the process: geometry, meshing, and setup. I won't go into too much details for the tedious parts but enough for you to figure it out. The most important part is the setup. Here I will go through each step thoroughly. 

I kind of used an unorthodox methodology for my simulations. I used ANSYS Workbench to create my geometry and meshing, but would open up Fluent directly rather than using Fluent through Workbench. After making my geometry and meshing, I would save the mesh file and import that to Fluent directly. This worked best for me because I could try a lot of different setups easily without having to save a lot of Workbench projects.

## Geometry
To create my channel geometry I used SpaceClaim through Workbench. Since we are dealing with membranes, we need work with micrometers to model the true thickness of the membrane. I used 40 micrometers for the membrane thickness but I should have used a smaller value that was cited in the paper. Nevertheless, you need to tell Workbench that you want to work at a smaller scale. You can do so by clicking the "Units" tab in the top left and clicking the last metric unit system that includes $\mu m$. 

![image](https://github.com/scanton1/CRO-Project/assets/168468147/96d153e5-1f3b-4687-8e30-72b339cfa7fe)

Now you can open up SpaceClaim or a different CAD software. Once you open up SpaceClaim, you need to click "File" in the top left corner and select "SpaceClaim Options". This will pop open the options menu and in here you should select "Units". From there you should select the dropdown menu for "Length" and "Micrometers" should be selected. At this point you can create your geometry with micrometers as your length unit. This inconvenient when creating the larger parts of your domain but it is the only method I know of to get the small scale needed for modeling a membrane.

The channel geometry I used has a feed height of 4.3 mm and length of 220 mm. Technically the feed height is bigger because I have three explicit zones, each 40 micrometers, that are part of the feed. This is insignificant and can be ignored honestly. These explicit zones are created for the UDF to work. The picture referenced in the UDF section shows the layout of the geometry. There is a membrane adjacent zone that is used to apply the UDF that assigns the species mass fraction, and there are two zones above, cell zone A and cell zone B, that are used for that UDF. Lastly, there is a zone for the membrane itself at the bottom of the geometry. The two pictures below show you the whole domain and a zoomed in view of the domain. 

![image](https://github.com/scanton1/CRO-Project/assets/168468147/b71ce9fc-9ee4-44fa-80a0-2f1c0970c232)

## Meshing
When you open up the meshing software, you must first make sure only the boundary face zones are outlined in red. If an interior face is outlined red then something went wrong when making your geometry. This is a result of how the overlapping zones are related in the topology. To fix this you need to go back to SpaceClaim and click on the parent node under the "Structure" menu on the left hand side. When this is selected there will be a section under properties called "Share Topology". In this drop down menu change it to "Share". This should fix the problem seen in the meshing software. 

When it comes to the meshing, I used a simple structured mesh composed of quad cells with inflation above the cell zone B and at the top wall. I make it so the explicit zones I made are only one cell thick. This ensures that the code I wrote does exactly what I want it to. In the future you could make the mesh more fine over those regions, but that would mess up the the discretization for the $C_{m}$ calculations. You could alter the code accordingly. When you create the mesh you need to make sure to create named selections for at least the following sections. It will likely change as the geometry becomes more complicated though.  

![image](https://github.com/scanton1/CRO-Project/assets/168468147/c26871b0-6afd-447b-a00a-39666686b2d1)

A section of the channel mesh can be seen below. Notice that there are four zones at the bottom all with the same thickness. Also you can see that the inlet is outlined in blue. A flaw with the simulation I've done is that the sides of the explicit zones above the membrane are considered as walls. Technically they should be part of the inlet or feed outlet but I considered them walls out of convenience. 

![image](https://github.com/scanton1/CRO-Project/assets/168468147/f2eda832-e6ff-494e-a979-18fc6a816018)

## Setup
Once satisfied with your mesh, you can open up Fluent and import your mesh.

![image](https://github.com/scanton1/CRO-Project/assets/168468147/c03abad9-0c87-492a-9eb3-2dba3a573333)

With Fluent open, it's time to setup the simulation. The first thing you can do is load in your UDFs. You can do this by going to the "User-Defined" menu along the top of Fluent and clicking "Functions". Clicking "Managed..." allows you to see existing libraries loaded in. When you make changes to your code you have to unload the library and recompile everytime. To compile and load in your UDF library, click "Compiled..." and use the GUI to add your source files and header files. For this code, there is only a single source file named "cro.c". Make sure the file is selected such that "Source Files [1/1]", then you can select "Use Built-In Compiler" and hit "Build". The compiler then tries to compile your code. Make sure to check the console for any errors while compiling. If there aren't any, you can click "Load". You should see in the console each UDF that was loaded in and ready for use.  

https://github.com/scanton1/CRO-Project/assets/168468147/e6162c97-44c1-4103-8f15-bea6c84ff3e9

With the UDFs loaded in, you can now use any of the functions and do the rest of the setup. We can first explore how to hook some of our UDFs and execute them. In the same "User-Defined" menu, there are the "Function Hooks..." and "Execute on Demand..." menus. The former can be used to hook/apply certain types of UDFs. The latter can be used to execute UDFs at any time. I attached a video below where I show you how to use these menus. 

https://github.com/scanton1/CRO-Project/assets/168468147/6057b0f4-2d7c-4e47-a633-7599ea4826a8

Let's go into the setup now. Under the "General Page", we select "Pressure-Based" and "Transient". 

![image](https://github.com/scanton1/CRO-Project/assets/168468147/fbde6280-ae06-4b0d-a3a2-9e803d064dc9)

Under the "Models" tab, select "Viscous" and choose "Laminar" for the channel problem. If you are doing simulations with turbulence, pick the appropriate model. Still under the "Models" tab, now select "Species" and turn on the "Species Transport" model. This will automatically default to the base mixture template.  

![image](https://github.com/scanton1/CRO-Project/assets/168468147/41ca7de2-f764-42c5-a62a-cab47778b7ff)

Now we need to create the salt water mixture. Double click the "Materials" tab and a page should pop up showing all the defined materials. On this page double click on "Fluid". A new page will open allowing you to define a new fluid. 

https://github.com/scanton1/CRO-Project/assets/168468147/19cf7cca-d813-40a6-ae69-b4f6adebd225

With both species individually created, we can now define the mixture to be salt water. 

https://github.com/scanton1/CRO-Project/assets/168468147/85472fd5-df42-442f-869f-65841aaef540

Under the "Cell Zone Conditions" menu there should be a list of fluid cell zones. We can alter how Fluent treats these zones by editing the cell zone. For my case, I apply a porous zone to the membrane zone and fixed values to the membrane and the membrane adjacent zones.  

https://github.com/scanton1/CRO-Project/assets/168468147/a33458ea-37c8-40ec-aada-69de7ac83b10

https://github.com/scanton1/CRO-Project/assets/168468147/8a5ef77c-94c9-40f9-b58f-0a483e2eebe8

https://github.com/scanton1/CRO-Project/assets/168468147/610e9821-2f6a-4b89-8e23-c1256ae571aa

Under the "Boundary Conditions" menu we can supply our boundary conditions. 

https://github.com/scanton1/CRO-Project/assets/168468147/5c604ec1-1031-49b1-9464-6badc67fad7a

https://github.com/scanton1/CRO-Project/assets/168468147/a16f0773-d82f-436b-a7d1-db0dc2487e31

Lastly, we can define our solution methods. I chose to follow exactly what the paper I referenced did. 

https://github.com/scanton1/CRO-Project/assets/168468147/d34d51ae-394b-4eba-b72a-1e63a19d7bfa

That should cover the setup for the simulations, now we can briefly look at the solution following a 4-min flow time simulation.

https://github.com/scanton1/CRO-Project/assets/168468147/0a73fa7f-62ba-4a2a-9f76-10a1285e97d9

https://github.com/scanton1/CRO-Project/assets/168468147/2c58d576-f62a-422f-b07d-f4c699066d3e

## Future Work
Hopefully this overview was sufficient to introduce you to the UDFs I wrote and the setup for the channel simulations. Going forward a couple of things need to be done to advance this project:
1. Update the equations used for solution properties as functions of concentration. This equations were created for the sodium-sulfate, water solution and aren't compatible with seawater simulations.
2. Use a source term to simulate the rotation of the system. Or figure out a different way to do it.

I'm sure more things will come up than these two, but this is a good starting point at least. 

























