# Thermodynamics PS1
Project Submission: 
Thremodynamics ASTR 507 - PS 1
Anastastasios (Andy) Tzanidakis - atzanida@uw.edu


Table of Content: 
[TOC]


## Open Source Code Citation
The particle-box collision I independently developed for this assignment is publicly available on GitHub: [github.com/AndyTza/thermo_507/ps1](https://github.com/AndyTza/thremo_507). 

For this project I utilize the 


## Problem 1 
Here I will describe all the steps I followed. 

### Algorithm Overview

The main component of the algorithm comes from the initialization of the `Particle` class. The `Particle` class will hold the particle positions, velocities, radius, and mass. Furthermore, the `Particle` class will have a draw method that will allow us to draw the particle in any Cartesian position: 

```python
class Particle():
    """
    Initialize a particle class with:
    position (2D), velocity (2D), radius and mass

    Input
    -----
    pos_x (float or int): Initial position in the x-axis direction
    pos_y (float or int): Initial position in the y-axis direction
    vx (float or int): Initial velocity in the x-axis direction
    vy (float or int): Initial velocity in the y-axis direction
    rad (float or int): Radius of each particle
    """
    def __init__(self, pos_x, pos_y, vx, vy, rad, mass):
        self.pos_x = pos_x
        self.pos_y = pos_y
        self.vx = vx
        self.vy = vy
        self.rad = rad
        self.mass = mass
    
    def draw_particle(self, window, color=(0, 0, 0)):
        """Draw a circular particle given a window and color.

        Args:
            window ([pygame.display]): Display class for window
            color (tuple, optional): Set the color of the particle. Defaults to a black color in RGB (0, 0, 0)

        Returns:
            [pygame.draw.circle]: Returns pygame circle object drawn
        """
        return pygame.draw.circle(window, color, (self.pos_x, self.pos_y), self.rad)
```

To draw the simulation for this problem I used the Python module `pygame` to draw a GUI background and add the particle shapes.`pygame` runs in a conditional while statement where the many functions of the particle collisions occur:

```python=
def simulate(N=100, stop_sim=50):
    """Initialize the pygame window & run simulation. 
    
       Input 
       -----
       N (int): Number of particles we want to simulate
       stop_sim (int): After how many iterations do you want the simulation to stop
    """
    
    # Begin pygame clock 
    clock = pygame.time.Clock()
    run = True
    
    # Initialize particles (see next sections)
   
    while run:
        clock.tick(FPS) # keep to the limited FPS
        draw_window_bkg() # draw window background -- need to refresh each 
    
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                run = False
                           
        # Run particle actions (i.e particle-wall or particle-particle collision)
        
        # Drawing & keeping track of particles
        k = 1
        for particle in particles_N:
            print (f"{k}-ID: x: {particle.pos_x}, y:{particle.pos_y}, Vx {particle.vx}, Vy:{particle.vy}")
            particle.draw_particle(WIN, color=(0, 0, 0))
            k+=1

        # Update pygame display!
        update_display()

        #print (f"Time interval: {delta_t} seconds")
        delta_t += 1 
        if delta_t>stop_sim:
            run = False

```

The background is controlled by these global parameters: 
- 800x800 pixel GUI box (highlighted in line 30 from the code above)
- 30 Frames per second simulation speed


In the following sections I will describe in more detail how the particles collide and move through the confined box.


### Particle Spawning & Initial Parameters


:::danger
**Important**: Avoiding particle-particle spawn issue. Occas
:::







### Handling a Particle-Wall Collision

The simplest scenario in an elastic collision is to consider the classical particle-wall collision. This can be achieved by examining the (x,y) position of each particle. We need to remember however that our particles have some radius ($r_p$). 

If the position in either axis plus the radius has passed the wall boundary condition then a collision has occurred, we will flip the sign of each velocity component: 


\begin{cases}
Collision \space \hat{x},& \text{if } 0 \leq x+r_{p}\geq X_{wall}\\
Collision \space \hat{y},& \text{if } 0 \leq y+r_{p}\geq Y_{wall}
\end{cases}

```python

def check_wall_collision(particle, axis='x', wall_boundary=800):
    """Will check if a particle-wall collision has occurred upon a specified axis.
       
       Input
       -----
       particle (class): Class particle
       axis (str): Axis of particle-wall collision (default to x-axis)       
       wall_boundary (float): Wall boundary (i.e from 0 to wall boundary)
        
       Output
       ------
       particle_wall_collision (bool): Returns boolean if collision has occurred
    """
    if axis=='x':
        if (particle.pos_x - particle.rad) <=0 or (particle.pos_x -particle.rad) >= wall_limit:
            return True
        else:
            return False
    if axis=='y':
        if (particle.pos_y - particle.rad) <=0 or (particle.pos_y -particle.rad) >= wall_limit:
            return True
        else:
            return False 
```

With the `check_wall_collision` function, now for each particle iteration, I check to see if a particle collision has occurred for both unit axes. If the collision occurs, the velocity get updated by the same magnitude but with a negative sign:

```python
# handle wall collisions
for ppp in particles_N:
    check_x = check_wall_collision(ppp, axis='x')
    check_y = check_wall_collision(ppp, axis='y')

    if check_x==True:
        ppp.vx = -ppp.vx # flip sign of x-component velocity
    if check_y==True:
        ppp.vy = -ppp.vy # flip sign of y-component velocity

```

### Handling the Particle-Particle Collision
```python

def particle_2_particle_collision(p1, p2):
    """Describe the particle to particle collsion"""
```





