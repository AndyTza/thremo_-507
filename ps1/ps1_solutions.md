# Thermodynamics PS1
Project Submission: 
Thremodynamics ASTR 507 - PS 1
Anastastasios (Andy) Tzanidakis - atzanida@uw.edu


Table of Content: 
[TOC]


## Open Source Code Citation
Here I present the ideal gas particle-box collision I independently developed in Python. The code is publicly available on GitHub: [github.com/AndyTza/thermo_507/ps1](https://github.com/AndyTza/thremo_507)


## Problem 1 - Ideal Gas Simulation
In this problem I will briefly describe how the algorithm works. In the animation below, you can see what a typical run looks like using this algorithm. 

|<img src="https://i.imgur.com/eGo830d.gif" alt="Image" width="300" height="300" style="display: block; margin: 0 auto" />|
|:--:| 
| Figure 1. -- Animation of ideal gas in a box and demo of the current code. In the following sections I will be describing in detail on how such calculations are carried out.|

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

```python
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

:::warning
**Important**: 

Throughout my simulations I use the following conditions: 
- 1000x1000 pixel GUI box (highlighted in line 30 from the code above)
- 30 frames per second (FPS) simulation speed
:::



### Particle Spawning & Initial Parameters

:::danger
**Important**:
Avoiding particle-particle spawn issue. In a finite box, occasionally two particles can spawn at roughly the same position. To avoid this issue, I initialize the position of all my (x,y) particle vectors in an evenly spaced grid such that the separation is always at least $\geq r_{p}$. After only a few iterations I found that this does **not** impact any of the distributions of the particle phase-space distributions.
:::

The particle positions are spawned by a 10x10 grid such that they keep a **minimum** separation of $r_{p}$ as seen by this exampe:

|<img src="https://i.imgur.com/476NH3C.png" alt="Image" width="300" height="300" style="display: block; margin: 0 auto" />|
|:--:| 
| Figure 2. -- (3x3) snapshot of my initial particle positions at time=0. We can see from thsi figure that the starting positions are at least $\geq r_{p}$. |

For the initial velocity ($v_0$), I assume a value of 200 px/s for both the xy components but with random sign. This enables all particles to have the same kinetic energy. After experimenting with a variety of initial velocities, I find that this value avoids the least the effects of tunnelling[^1]. 

For particle mass I assume that all particles will have the same mass (m=10). Finally, all particles will also have the same radius ($r_{p}$) that is approximately 2% of the container size. 


[^1]: When the velocity of a particle is too fast (i.e faster than the average processing time per step) you can often encounter unwanted effects such as particles not colliding. 






### Handling a Particle-Wall Collision

The simplest scenario in an elastic collision is to consider the classical particle-wall collision. This can be achieved by examining the (x,y) position of each particle. 

If the position in either axis plus the radius has passed the wall boundary condition then a collision has occurred, we will flip the sign of each velocity component: 


\begin{cases}
Collision \space \hat{x},& \text{if } 0 \leq x+r_{p}\geq X_{wall}\\
Collision \space \hat{y},& \text{if } 0 \leq y+r_{p}\geq Y_{wall}
\end{cases}

```python
def check_wall_collision(pt, axis='x', wall_limit=width):
    """Check for wall collision. Returns if a wall collision has occured on a specified axis. Returns bool if collision has occured.

        Input
        -----
        pt (class; particle): Particle class
        axis (str): Axis of collision, currently supports 'x' and 'y' collisions
        wall_limit (float): Wall limit assuming a box window

        Output
        ------
        bool; Will return the bool of the collision status

    """
    if axis=='x':
        if pt.pos_x-pt.rad <= 0 or pt.pos_x+pt.rad >= wall_limit:
            return True
        else:
            return False
    if axis=='y':
        if pt.pos_y-pt.rad <= 0 or pt.pos_y+pt.rad >= wall_limit:
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
In the code above, `particles_N` are a list of particle classes. For each time step, I check for wall collisions. If a wall collision has occurred then we will reflect the velocity in each component but with opposite sign.  


### Handling the Particle-Particle Collision
I also take into consideration the particle-particle collision. In an elastic collision scenario we know that the momentum and kinetic energy will be **conserved**. 

First, I check if any unique particle pair has collided. To do this, I calculate the Euclidean distance for each unique pair and declare the status of the collision if: 

\begin{equation}
d_{2,1} \leq 2r_p - C
\end{equation}

and 

\begin{equation}
<\vec{v_1}-\vec{v_2}, \vec{r_{1}}-\vec{r_2}> \space \lt  0 
\end{equation}

where C is some boundary constant for adjusting the spacing between particles.

[^2]: A nearest neighbor algorithm can help speed this process up.


After a particle-particle collision has been detected I apply an elastic velocity exchange. For the particle-particle collision, I chose to express the collision in an angle-free format[^3]:

[^3]: [Angle-free elastic wiki](https://wikimedia.org/api/rest_v1/media/math/render/svg/14d5feb68844edae9e31c9cb4a2197ee922e409c)

\begin{equation}
\vec{v'_{1}} = u_{1,x} - \frac{2m_2}{m1+m2} \frac{<v_{1,x}-v_{2,x}, r_{1}-r_{2}>}{||x_1-x_2||^2} (r_1-r_2) \\
\vec{v'_{2}} = u_{2,x} - \frac{2m_1}{m1+m2} \frac{<v_{2,x}-v_{1,x}, x_2-x_1>}{||x_2-x_1||^2} (x_2-x_1)
\end{equation}

where $\vec{v'_{1}}$, $\vec{v'_{2}}$ are the updated velocities after collision. 

In practice, this is computed with each unique pair of particles: 

```python=
# handle particle-particle collisions
_off = 0
for i in range(0, len(particles_N)-1):
    jj = np.arange(_off+1, len(particles_N))
    for j in jj:
        
        # Check for collisions
        _cp = check_particle_particle_collision(particles_N[i], particles_N[j], boundary=3) #

        pt_1 = particles_N[j] # i-th particle
        pt_2 = particles_N[i] # j-th particle

        # positional vector
        vec_1 = np.array([pt_1.pos_x, pt_1.pos_y])
        vec_2 = np.array([pt_2.pos_x, pt_2.pos_y])

        # velocity vector
        vel_1 = np.array([pt_1.vx, pt_1.vy])
        vel_2 = np.array([pt_2.vx, pt_2.vy])

        m1, m2 = pt_1.mass, pt_2.mass
        M = m1 + m2

        # normalized distance!
        cart_dist = ((vec_1 - vec_2)[0]**2 + (vec_1 - vec_2)[1]**2)

        if _cp==True and np.dot(vec_1 - vec_2, vel_1 - vel_2)<0:

            # Collision updated velocity!
            vel1_new = vel_1 - (2*m2)/M * np.dot(vel_1 - vel_2, vec_1 - vec_2) * (vec_1 - vec_2)/cart_dist
            vel2_new = vel_2 - (2*m1)/M * np.dot(vel_2 - vel_1, vec_2 - vec_1) * (vec_2 - vec_1)/cart_dist

            pt_1.vx = vel1_new[0]
            pt_1.vy = vel1_new[1]
            pt_2.vx = vel2_new[0]
            pt_2.vy = vel2_new[1]

            # Now apply equations of motion to the actual positions (i.e make them move!)
            pt_1.pos_x = pt_1.pos_x + pt_1.vx*(1/FPS)
            pt_1.pos_y = pt_1.pos_y + pt_1.vy*(1/FPS)

            pt_2.pos_x = pt_2.pos_x + pt_2.vx*(1/FPS)
            pt_2.pos_y = pt_2.pos_y + pt_2.vy*(1/FPS)

        else:
            continue
    _off += 1

```


### Initial Versus Final Velocity Distributions

Here I test the initial and final velocity distribution for my ideal gas simulator. I begin with initial conditions: 
- Wall boundary: 1000$px$ x 1000$px$ periodic boundaries
- Initial velocity: ||$u_{0x}$, $u_{0y}$|| = 200 $px$/s (with random sign)
- Particle mass: m=10
- End simulation at $\delta t$=10$^{4}$ seconds


|<img src="https://i.imgur.com/l1KrNaF.jpg" alt="Image" style="display: block; margin: 0 auto" />|
|:--:| 
| Figure 3. -- Initial versus final position (left panels) and velocity distribution (right panels). The top row shows the initial position with the velocity vector component highlighted with its constant distribution in initial velocity (top right). The bottom panels show the final condition of the simulator after $10^{4}$ seconds. We see that the positions of the particles are randomly dispersed while their velocity distribution has a characteristic Maxwellian distribution.|

## Problem 2 - Maxwellian

### Derivation of the Maxwellian Distribution

Let us begin by assuming the functional form of the Maxwellian distribution in 2D: 

\begin{equation}
f(E) = \alpha e^{-\frac{E}{k_{B}T}}
\end{equation}

#### Solving the $\frac{dN}{dv}$ Distribution: 

Generally, we know that:
\begin{equation}
\frac{dN}{dv}= f(E) = \alpha e^{-\frac{E}{k_{B}T}}
\end{equation}

\begin{equation}
dN = f(E) dv
\end{equation}

For our 2D particle box simulation, we an assume that each particle will have a velocity and angle ($\phi$) component. In the projected space this can convert $d\vec{v}$=$v$d$u$d$\phi$:

\begin{equation}
dN = f(E) v dv d\phi
\end{equation}

By normalizing the area under the curve we can take:

\begin{equation}
f(E)v dv d\phi =1
\end{equation}

Now, let's expand and find parameter $\alpha$:


\begin{equation}
\int_0^{\infty} \int_0^{2\pi} \alpha e^{-\frac{E}{k_{B}T}} v dv d\phi = 1
\end{equation}

Generally, for the $n_{i}^{th}$ particle will have energy $\frac{1}{2} m v^2$:

\begin{equation}
\int_0^{\infty} \int_0^{2\pi} \alpha e^{-\frac{E}{k_{B}T}} v dv d\phi = 1
\end{equation}

\begin{equation}
\int_0^{\infty} \alpha e^{-\frac{E}{k_{B}T}} v dv \int_0^{2\pi} d\phi = 1
\end{equation}


\begin{equation}
\int_0^{\infty} \alpha e^{-\frac{E}{k_{B}T}} v dv 2\pi = 1
\end{equation}

First, let's expand our energy term and re write:

\begin{equation}
\int_0^{\infty} \alpha v 2\pi exp(-\frac{mv^2}{k_B T}) = 1
\end{equation}

\begin{equation}
\alpha = \frac{1}{2 \pi \int_0^{\infty} v\space exp(\frac{-mv^2}{2k_{B}T})dv}
\end{equation}

Now let's solve the following integral

\begin{equation}
\int_0^{\infty} v\space e^{\beta v^2}dv \sim -\frac{1}{2 \beta^2}
\end{equation}

where $\beta=\frac{-m}{k_{B}T}$

By pluggning this back into the $\alpha$ equation, we find:

\begin{equation}
\alpha = \frac{m}{2\pi k_{B} T}
\end{equation}

Finally, we can re-plug in our $\alpha$ term to our original equation to get:

\begin{equation}
\frac{dN}{dv} = \bigg{(}\frac{m}{2\pi k_{B} T} \bigg{)}e^{-\frac{E}{k_{B}T}}
\end{equation}

#### Solving the $\frac{dN}{dE}$ Distribution: 

Similarly to my derivation above, now we can find the $\frac{dN}{dE}$ distribution function. Let us begin with:

\begin{equation}
\frac{dN}{dE} = \gamma e^{-\frac{E}{k_{B}T}}
\end{equation}

we can rearrange our equation:

\begin{equation}
dN = \gamma e^{-\frac{E}{k_{B}T}} dE
\end{equation}

Note that we can express the change in energy and write it in terms of velocity and impact angle per particle $dE=mvdv=mvvdvd\phi=mv^2dvd\phi$:

\begin{equation}
dN = \int_0^{\infty}\int_0^{2\pi} \gamma e^{-\frac{E}{k_{B}T}} mv^v dv d\phi
\end{equation}

I will now normalize the integral and solve for $\gamma$:

\begin{equation}
\int_0^{\infty}\int_0^{2\pi} \gamma e^{-\frac{E}{k_{B}T}} mv^2 dv d\phi = 1
\end{equation}


\begin{equation}
2\pi \gamma m \int_0^{\infty}  v^2 e^{-\frac{mv^2}{2k_{B}T}} dv d\phi = 1
\end{equation}

Let's now solve the following integral assuming $\delta=-\frac{m}{2k_{B}T}$:

\begin{equation}
\int_0^{\infty} v^2 e^{\delta v^2} dv = -\frac{\sqrt{\pi}}{4\delta^{3/2}}
\end{equation}

After plugging everything back into our equation and solving for $\gamma$:

\begin{equation}
\gamma = \frac{1}{\sqrt{\frac{2}{m}} (\frac{k_{B}T}{\pi})^{\frac{3}{2}}}
\end{equation}

Finally, we can take our $\gamma$ factor and put it back to our original equation: 

\begin{equation}
\frac{dN}{dE} \sim \frac{1}{\sqrt{\frac{2}{m}} (\frac{k_{B}T}{\pi})^{\frac{3}{2}}} exp(-\frac{mv^2}{2k_{B}T})
\end{equation}

### Distribution in Velocity - Maxwellian vs. Ideal Gas Particle-Box

Using the equations from above, we can now plot the expected Maxwellain distribution and see how well it overall agrees with our observed distribution. 

In the following figure, I show the theoetical (gray solid) fitted Maxwellian cumulative density function and observed data. The left figure shows the distribution for one simulation, and the right for 10 simulations ending at $10^{4}$ seconds. **I found strong agreement between the classical Maxwellian distribution and the observed velocity dispersion.**

|<img src="https://i.imgur.com/7KkxAUX.jpg" alt="Image"  style="display: block; margin: 0 auto" />|
|:--:| 
| Figure 4. -- Comparison of the cumulative density function of a theoretically derived Maxwellain distribution versus the velocity distribution observed with the ideal gas simulator. The left panel shows the CDF on a single simulation. Similarly, the right figure shows the CDF distribution for ten stacked simulations.|

### Changing Mass Parameter 

For the following analysis, I split the particle mass. Half of the particles were assigned with mass $m=10$ and the other half with 10$*m$. For the simulation, I ran the simulation 50 times with an ending time of 2x$10^3$ seconds (assuming that most particles have reached equilibrium). After stacking all velocities from each simulation, I fit a Maxwellian on the "lighter" particles (solid line; Fig. 4) and "heavier" particles (dashed line; Fig. 4). 

**Overall, I found that the velocity distribution now was primarily dominated by two Maxwellian components (i.e the heavier mass and lighter mass components).** For the estimated energy distribution (i.e $\frac{1}{2}mv^2$) I found the distribution to have a very sharp rise and a longer tail at higher energies (likely caused by the more massive particles).

|<img src="https://i.imgur.com/DBOAkVR.jpg" alt="Image"  style="display: block; margin: 0 auto" />|
|:--:| 
| Figure 5. -- Probability density function of the velocity distribution for 50 stacked simulations. Half of the particles in this simulation are $10m$ the original mass particle. We find that the distribution can be well described by a mixture of two Maxwellians at two different velocity components.|


## Problem 3 - Relaxation Timescales

### Investigating Possible Relaxation Timescale Metrics

For T time steps in our simulation, we can track the dispersion in velocity for each time step. At some relaxation time ($T_{R}$), we expect that the average dispersion in velocities to be roughly Gaussian and be perturbed by a central dispersion of the system. 

In the figure below I show the standard deviation of the velocity for each time step in a simulation. I track $\sigma_{v}$ for at least 10,000 sec. 



|<img src="https://i.imgur.com/SerBmbM.jpg" alt="Image"  style="display: block; margin: 0 auto" />|
|:--:| 
| Figure 6. -- Time series of the standard deviation in velocity ($\sigma_{v}$) for a single simulation. We see that after some time $T_{R}$ we expect that our ideal gas simulation to reach an equilibrium where the $\sigma_{v}$ should average around a maximum value. We can use the stationarity of such time series feature to asses when our system has reached equilibrium. |

One possible metric to evaluate at each time is to estimate if cumulatively our time series has reached an equilibrium. I carry out an **Augmented Dickey–Fuller** (ADF) test[^2] to reject or accept the null hypothesis if the time series is **stationary** or **non-stationary**. For each time frame I estimate the cumulative ADF alongside the associated p-value. After t=1, I evaluate the median p-value of each ADF test. Once the velocity dispersion begins to reach it's relaxation time (i.e equilibrium), I found that overall the average p-value was $\leq$ 0.05 indicating that we can reject the null hypothesis that the time series does **not** have a unit root (meaning it's stationary).

[^2]: The ADF test essentially tests the time series data to reject or accept the null hypothesis if there is a unit root. In simple, a unit root will dictate if the time series has reached a stationary condition that is near the expected value (i.e equilibrium).

For this test of 1000 time steps, I found that this particular system reached a stationary status around 500 sec. 

|<img src="https://i.imgur.com/6ATBSeG.jpg" alt="Image"  style="display: block; margin: 0 auto" />|
|:--:| 
| Figure 7. -- Tracking the p-value for the ADF test for each time step of our simulation. We find that after some time the p-values reach a threshold that is p$\leq$ 0.05 indicating time time at where our time series velocity dispersion has reached its equilibrium. |


### Relaxation Timescale per Particle
We want to derive the relaxation timescale per particle with the following parameters: 

\begin{equation}
u_0, n, \sigma
\end{equation}

where $u_{0}$ ($ms^{-1}$) is the initial particle velocity, $n$ ($m^{-2}$) the number density, and cross section of particle $\sigma$ ($m$). If we apply dimensional analysis we can find the relaxation timescale: 

\begin{equation}
T_{R} \sim \frac{1}{u_0 n \sigma }
\end{equation}

### Metric Relaxation Timescale versus Derived Relaxation Timescale

In this last question, I compare the relaxation timescale of the derived equation versus the ADF test metric. To vary the cross section of each particle I change the particle size by a factor of 5$r_{p}$. For 50 stacked simulations, I found a somewhat good agreement with the the derived metric. I suspect that the first few simulated runs do not match quite well the derived relations because the dispersion in velocity has a much higher variation in velocity even when it reached its equilibrium state. This will cause the ADF test to incorrectly identify the feature as stationary. 

|<img src="https://i.imgur.com/rF1dCAu.jpg" alt="Image"  style="display: block; margin: 0 auto" />|
|:--:| 
| Figure 6. -- Expected timescale relaxation time (derived from the relaxation time above) versus our invented relaxation time ADF metric as a function of particle cross section.|

