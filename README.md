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

### Content
* `showcase_poisson_bolzmann.ipynb`: A working example
* `generate_structure.py`: Sampling and plotting
* `poisson_bolzmann_distribution.py`: generate potential and densities

### Usage
If you can't get the notebooks to run of the box, set up an environment for them:
```bash
mkdir env
virtualenv env
source env/bin/activate
ipython kernel install --user --name=new_kernel
jupyter notebook showcase_poisson_bolzman_distribution.ipynb
```
And then choose `new_kernel` in the top right dropdown menu as a kernel.
