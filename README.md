# Generate atomic distribution
Generate an MD input file containing a number of atoms in solution, which approximate a given distribution in space.

### Background
In order to understand lubrication better, we simulate thin layers of lubricant on a metallic surface, solvated in water.
Different structures of lubricant films are created by varying parameters like their concentration and the charge of the surface.
The lubricant is somewhat solvable in water, thus parts of the film will diffuse into the bulk water.

![pic](https://i.ibb.co/Yh8DxVM/showpicture.png)

### Physical solution
Lubricant molecules are charged, and their distribution is described by the Poisson-Bolzmann equation.

### Simulating the structure
To incorporate these results into our lubricant simulation, we need to sample the correct distribution.
Then we need to transform these samples into the correct file formats and feed them into a Molecular Dynamics simulation.
