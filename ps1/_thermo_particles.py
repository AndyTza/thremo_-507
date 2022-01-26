from shutil import move
import numpy as np
import pygame 
import matplotlib.pyplot as plt
import os
from tqdm import tqdm
import random

import matplotlib.pyplot as plt
pygame.font.init()
time_stamp = pygame.font.SysFont('Times', 40)

# Pygame set-up
global width, height, WIN, CLR, PRT_CLR, FPS

width, height = 1000, 1000 
WIN = pygame.display.set_mode((width, height)) 
CLR = (255, 255, 255)
PRT_CLR = (0, 1, 0)
FPS = 60


class Particle():
    """
    Initialize particle position and velocities. 

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
        return pygame.draw.circle(window, color, (self.pos_x, self.pos_y), self.rad, width=0)

def draw_window_bkg(color=CLR):
    """Draw background according to the CLR argument"""
    WIN.fill(color)

def update_display():
    """Update display"""
    pygame.display.update()

def cart_distance(p1, p2):
    """Calculate the eucledian distance between particle 1 and particle 2."""   
    return np.sqrt((p2.pos_x - p1.pos_x)**2 + (p2.pos_y - p1.pos_y)**2)

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

def check_particle_particle_collision(particle_1, particle_2, boundary=100):
    """Given a pair of particles (of equal radius?!) check to see if they have collided. Retuns bool"""

    d = cart_distance(particle_1, particle_2) # cartesian distance in 2D
    
    if d<=(particle_1.rad + particle_2.rad - boundary):
        return True
    else:
        return False


def spawn_particles(N0, vmin=-100, vmax=100, rad=0.011*height, m0=10, seed=42):
    """ This function will spawn N particles in WIN box. It will initialize the position, velocity, mass and size of each particle.

    Args:
        N0 (int): Number of spawned particles
        particle_sigma (float): Dispersion in particle position
    """

    _step_rec = np.sqrt((width**2)/N0) #reccomended step to evenly space particles

    X, Y, offset = [], [], 50 # apply some vertical and horizontal offset...
    for i in np.arange(0, width, step=_step_rec):
        for j in np.arange(0, height, step=_step_rec):
            X.append(i)
            Y.append(j)
    pos_rand_x, pos_rand_y = np.array(X), np.array(Y)
    pos_rand_x, pos_rand_y = pos_rand_x + offset, pos_rand_y + offset
    #pos_rand_x, pos_rand_y = np.random.randint(5, width-5, N0), np.random.randint(5, width-5, N0)

    np.random.seed(seed)
    all_particles = [Particle(i, j, 6*random.choice((-1, 1)), 6*random.choice((-1, 1)), rad, m0) for i, j in zip(pos_rand_x, pos_rand_y)]

    return all_particles

def simulate(N=10, stop_sim=50, **kwargs):

    clock = pygame.time.Clock()
    run = True
    
    # Initialize particles
    particles_N = spawn_particles(N0=N, **kwargs)
    delta_t = 0 

    while run:
        clock.tick(FPS) # keep to the limited FPS
        draw_window_bkg() # draw window background -- need to refresh each 
    
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                run = False

        # Update velocity
        for pp in particles_N:
            pp.pos_x = pp.pos_x + pp.vx*(1/FPS)
            pp.pos_y = pp.pos_y + pp.vy*(1/FPS)

        for ppp in particles_N:
            check_x = check_wall_collision(ppp, axis='x')
            check_y = check_wall_collision(ppp, axis='y')

            if check_x==True:
                ppp.vx = -ppp.vx
                ppp.pos_x = ppp.pos_x + ppp.vx*(1/FPS)
            if check_y==True:
                ppp.vy = -ppp.vy
                ppp.pos_y = ppp.pos_y + ppp.vy*(1/FPS)


        # handle particle-particle collisions
        _off = 0
        for i in range(0, len(particles_N)-1):
            jj = np.arange(_off+1, len(particles_N))
            for j in jj:

                _cp = check_particle_particle_collision(particles_N[i], particles_N[j], boundary=3) # check if collision has occured
                
                # If particle collision has occured... perform a velocity change
                pt_1 = particles_N[j] # i-th particle
                pt_2 = particles_N[i] # j-th particle

                m1m2 = pt_1.mass + pt_2.mass # mass combined term
                
                # x-hat component
                x1x2 = pt_1.pos_x - pt_2.pos_x
                x2x1 = pt_2.pos_x - pt_1.pos_x
                v1v2 = pt_1.vx - pt_2.vx
                v2v1 = pt_2.vx - pt_1.vx

                y1y2 = pt_1.pos_y - pt_2.pos_y
                y2y1 = pt_2.pos_y - pt_1.pos_y
                v1v2y = pt_1.vy - pt_2.vy
                v2v1y = pt_2.vy - pt_1.vy

                if _cp==True and x1x2!=0 and x2x1!=0 and y1y2!=0 and y2y1!=0:

                    ## THE ISSUE IS HERE!!!


                    pt_1.vx = pt_1.vx - ((2*pt_2.mass)/(m1m2)) * ((np.inner(v1v2,x1x2))/np.linalg.norm(x1x2)**2) * (x1x2)
                    pt_2.vx = pt_2.vx - ((2*pt_1.mass)/(m1m2)) * ((np.inner(v2v1,x2x1))/np.linalg.norm(x2x1)**2) * (x2x1)

                    pt_1.vy = pt_1.vy - ((2*pt_2.mass)/(m1m2)) * ((np.inner(v1v2y,y1y2))/np.linalg.norm(y1y2)**2) * (y1y2)
                    pt_2.vy = pt_2.vy - ((2*pt_1.mass)/(m1m2)) * ((np.inner(v2v1y,y2y1))/np.linalg.norm(y2y1)**2) * (y2y1)

                    # Now apply equations of motion to the actual positions (i.e make them move!)
                    pt_1.pos_x = pt_1.pos_x + pt_1.vx*(1/FPS)
                    pt_1.pos_y = pt_1.pos_y + pt_1.vy*(1/FPS)

                    pt_2.pos_x = pt_2.pos_x + pt_2.vx*(1/FPS)
                    pt_2.pos_y = pt_2.pos_y + pt_2.vy*(1/FPS)

                else:
                    continue
            _off += 1

        k = 1
        for particle in particles_N:
            print (f"{k}-ID: x: {particle.pos_x}, y:{particle.pos_y}, Vx:{particle.vx}, Vy:{particle.vy}")
            particle.draw_particle(WIN, color=(0, 0, 0))
            k+=1

        # Let's draw a rectangle at the edge of the window display
        pygame.draw.rect(WIN, (0, 0, 0), (0, 0, 200, 50))
        _text = time_stamp.render(f"Time: {delta_t}", 1, (255, 0, 0	))
        WIN.blit(_text, (1, 1))
        
        # Update display!
        update_display()

        # update everything
        #for particle in particles_N:
            #particle.pos_x = particle.pos_x
            #particle.pos_y = particle.pos_y
            #particle.vx = particle.vx
            #particle.vy = particle.vy

        #print (f"Time interval score: {delta_t} seconds")
        delta_t += 1 
        if delta_t>stop_sim:
            run = False
            return particles_N


# Run the simulation here!
if __name__=="__main__":
    ss_1 = simulate(N=100 ,stop_sim=5000, vmin=-150, vmax=150)
    plt.figure(figsize=(5,5))
    for p in ss_1:
        plt.scatter(p.pos_x, p.pos_y, s=50, edgecolor='k')
        plt.quiver(p.pos_x, p.pos_y, p.vx, p.vy)
    for _ in range(2):
        plt.axvline(0, lw=5, color='tomato')
        plt.axhline(0, lw=5, color='tomato')
        plt.axvline(1000, lw=5, color='tomato')
        plt.axhline(1000, lw=5, color='tomato')
    plt.xlabel("X-pos", fontsize=18)
    plt.ylabel("Y-pos", fontsize=18)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlim(0, 1000)
    plt.ylim(0, 1000)
    plt.title("$\delta$t = 10$^{4}$ sec", fontsize=18)
    plt.savefig('figs/pos_final.jpg', format='jpg', dpi=400, bbox_inches='tight')
    plt.show()

    np.random.seed(42)
    #v0x, v0y = np.random.randint(-150, 150, 100), np.random.randint(-150, 150, 100)
    #V0 = np.sqrt(v0x**2 + v0y**2)
    v = []
    for p in ss_1:
        v.append(np.sqrt(p.vx**2 + p.vy**2))
    v = np.array(v)
    
    plt.figure(figsize=(10,5))
    _ = plt.hist(v, alpha=0.5, color='tomato', edgecolor='red', lw=3, density=False, histtype='stepfilled', label='Final')
   #_ = plt.hist(V0, alpha=0.5, color='c', edgecolor='k', bins=np.arange(0, 400, step=25), ls='--', density=True, histtype='stepfilled', label='Inidial')
    plt.ylabel("N", fontsize=20)
    plt.xlabel("$v$ [px/s]", fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.xticks(fontsize=18)
    plt.legend(fontsize=18)

    plt.title("$\delta$t = 10$^{4}$ sec", fontsize=18)
    plt.savefig('figs/dist_final.jpg', format='jpg', dpi=400, bbox_inches='tight')
    plt.show()
    
    

