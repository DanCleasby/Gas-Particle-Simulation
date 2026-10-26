import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from itertools import combinations


class Particles():
    
    
    def __init__(self, particle_number, radius, box_size):
        '''
        Initialisation of particle attributes and particle id pairs

        Parameters
        ----------
        particle_number : int
            number of particles
        radius : float
            radius of particles
        box_size : int
            box dimension

        Returns
        -------
        None
        '''
        
    
        self.radius = radius
        n = particle_number

        self.positions = np.random.uniform(self.radius, box_size - self.radius, (2, n))
        self.velocities = np.random.uniform(-500, 500, (2, n))
        
        self.ids = np.arange(n)
        
        ''' 2d array of id combinations'''
        self.id_pairs = np.asarray(list(combinations(self.ids,2)))
        
        
    def update_positions(self, dt):
        '''
        Euler method for updating position

        Parameters
        ----------
        dt : float
            timestep 

        Returns
        -------
        None
        '''
        
        self.positions = self.positions + self.velocities * dt
    
    
    
    
###############################################################################    
    
    
    
    
class Simulation():
    
    
    def __init__(self, particle_number, tfin, time_interval, radius, box_size):
        '''
        Stores main simulation attributes

        Parameters
        ----------
        particle_number : int
            number of particles
        tfin : float
            length of time for the simulation
        time_interval : float
            interval of each time step

        Returns
        -------
        None
        '''
        
        self.n = particle_number
        self.t = tfin
        self.dt = time_interval
        self.timesteps = round(self.t / self.dt)
        self.radius = radius
        self.box_size = box_size
        
      
    def sort_particle_ids(self, particles):
        '''
        lists combinations of all x and y coordinates, finds difference and then flattens to 1D array

        Parameters
        ----------
        particles : Particle class object

        Returns
        -------
        ndarray
            2d array of the particle ids which collide with eachother
        '''
        
        x_pairs = np.asarray(list(combinations(particles.positions[0], 2)))
        dx_pairs = np.diff(x_pairs, axis = 1).ravel() 
        
        y_pairs = np.asarray(list(combinations(particles.positions[1], 2)))
        dy_pairs = np.diff(y_pairs, axis = 1).ravel()
        
        d_pairs = np.sqrt(dx_pairs**2 + dy_pairs**2) 
        
        ''' finds the pairs where radius overlaps'''
        ids_collision = particles.id_pairs[d_pairs < 2 * particles.radius]
        
        return ids_collision
        
    
    def check_coll(self, particles):
        '''
        elastic collisions formula
        
        uses particle ids to check and compute collisions (reduces number of calculations)
    
        Parameters
        ----------
        particles : Particle class object
            
        Returns
        -------
        None
        '''
        
        ids_collision = self.sort_particle_ids(particles)
        
        ''' mask for particles which collide'''
        v1 = particles.velocities[:, ids_collision[:,0]]
        v2 = particles.velocities[:, ids_collision[:,1]]
        r1 = particles.positions[:, ids_collision[:,0]]
        r2 = particles.positions[:, ids_collision[:,1]]
    
        ''' np.sum instead of np.dot as easier with matrices '''
        dot1 = np.sum((v1-v2)*(r1-r2), axis=0)
        dot2 = np.sum((v2-v1)*(r2-r1), axis=0)
    
        ''' mask to check if particles are moving towards eachother, else particles could get stuck '''
        towards1 = dot1 < 0
        towards2 = dot2 < 0
        
        v1_final = v1 - dot1 * (r1-r2) / np.sum((r1-r2)**2, axis = 0)
        v2_final = v2 - dot2 * (r2-r1) / np.sum((r2-r1)**2, axis = 0)
        
        particles.velocities[:, ids_collision[towards1,0]] = v1_final[:, towards1]
        particles.velocities[:, ids_collision[towards2,1]] = v2_final[:, towards2]
        
        
    def check_reflection(self, particles):
        '''
        checks and computes collision with the box boundary

        Parameters
        ----------
        particles : Particle class object

        Returns
        -------
        None
        '''

        
        x_right_bool = particles.positions[0, :] > self.box_size - particles.radius 
        x_left_bool = particles.positions[0, :] < particles.radius
        y_up_bool = particles.positions[1, :] > self.box_size - particles.radius
        y_down_bool = particles.positions[1, :] < particles.radius
        
        positive_vx = particles.velocities[0, :] > 0
        negative_vx = particles.velocities[0, :] < 0
        positive_vy = particles.velocities[1, :] > 0
        negative_vy = particles.velocities[1, :] < 0
        
        ''' combining masks '''
        x_right_mask = x_right_bool & positive_vx
        x_left_mask = x_left_bool & negative_vx
        y_up_mask = y_up_bool & positive_vy
        y_down_mask = y_down_bool & negative_vy
            
        
        particles.velocities[0, x_right_mask] *= -1
        particles.velocities[0, x_left_mask] *= -1
        particles.velocities[1, y_up_mask] *= -1
        particles.velocities[1, y_down_mask] *= -1
        

        
    def simulate(self):
        '''
        runs the simulation using all other functions and saves data to a file

        Returns
        -------
        pos_array : ndarray
            3d array of form (timesteps, coordinates, particles)
        vel_array : ndarray
            3d array of form (timesteps, velocity, particles)
        '''
        
        particles = Particles(self.n, self.radius, self.box_size)
        
        pos_array = np.zeros([self.timesteps, 2, self.n])
        vel_array = np.zeros([self.timesteps, 2, self.n])
            
        pos_array[0] = particles.positions.copy()
        vel_array[0] = particles.velocities.copy()
        
        ''' checking stops particles getting stuck initially '''
        self.check_reflection(particles)
        self.check_coll(particles)
        
        for i in range(1, self.timesteps):
            
            particles.update_positions(self.dt)
            
            pos_array[i] = particles.positions.copy()
            vel_array[i] = particles.velocities.copy()
                
            self.check_reflection(particles)
            self.check_coll(particles)
            
        np.savez('MB_data_1000', pos=pos_array, vel=vel_array,n=self.n,
                 r=self.radius, box_size=self.box_size)
        
        return pos_array, vel_array
    
       
        
    
###############################################################################  
    
    
    
        
class Animation():
    
    def __init__(self, pos_array, vel_array, n, radius, box_size):
        '''
        takes data gathered from the Simulation class to animate the simulation and relevant graphs

        Parameters
        ----------
        pos_array : ndarray
            3d array of the form (timestep, coordinates, particles)
        vel_array : ndarray
            3d array of the form (timestep, velocity, particles)
        n : int
            number of particles
        radius : float
            radius of particles
        box_size : int
            dimensions of box

        Returns
        -------
        None
        '''
         
        self.pos_array = pos_array
        self.vel_array = vel_array
        self.n = n
        self.radius = radius
        self.box_size = box_size
        self.vel_mag_array = np.linalg.norm(self.vel_array, axis = 1)
    
    def animate_simulation(self):
        '''
        uses funcanimation and circle patches to animate the data set

        Returns
        -------
        None
        '''
        
        fig, ax = plt.subplots()
        
        ax.axis("scaled")
        ax.set_xlim(0, self.box_size)
        ax.set_ylim(0, self.box_size)
        ax.set_title('Balls')
        ax.set_xlabel('x axis')
        ax.set_ylabel('y axis', rotation = 0)
        ax.set_facecolor('black')
        
        ''' Normalize to use cmap to change color of particles based on velocity'''
        self.norm = Normalize(vmin=self.vel_mag_array.min(), vmax=self.vel_mag_array.max())
        self.cmap = cm.get_cmap('jet')
        
        self.circles = []
        
        for i in range(self.n):
            
            circle = plt.Circle((self.pos_array[i,0,0], self.pos_array[i,1,0])
                                , self.radius)
            ax.add_patch(circle)
            
            self.circles.append(circle)

        self.anim = FuncAnimation(fig, self.update_plot, frames=self.pos_array.shape[0],
                                  repeat=False, interval=5, blit=True)
                
        plt.show()
        
        
    def update_plot(self, i):
        """
        Update the plot for FunAnimation. Checks if velocity has changed, if so changes particle color.

        Parameters
        ----------
        i : int
            current frame number

        Returns
        -------
        list
            list of circle objects with updated centers

        """
        
        for j, circle in enumerate(self.circles):
            circle.set_center((self.pos_array[i,0,j], self.pos_array[i,1,j]))
            
            if self.vel_mag_array[i-1,j] != self.vel_mag_array[i,j]:
                
                ''' changes color based on current velocity'''
                color = self.cmap(self.norm(self.vel_mag_array[i,j]))
                circle.set_color(color)
        
        return self.circles
    
    def animate_graph(self):
        '''
        animates histogram of velocity data to show how the velocity distribution changes,
        this distribution eventually tends to the Maxwell-Boltzmann distribution

        Returns
        -------
        None
        '''
        
        self.max_velocity = self.vel_mag_array[:,:].max()
        self.bins = np.linspace(0, self.max_velocity + 10, 40)
        
        fig, self.ax = plt.subplots()
        
        self.ax.set_title('Velocity Distribution')
        self.ax.set_xlabel('Particle Velocity [m/s]')
        self.ax.set_ylabel('Number of Particles')
        self.ax.set_xlim(0, self.max_velocity + 10)
        self.ax.set_ylim(0, 120)
        
        self.y_ticks = np.arange(0,120, 10)
        self.ax.set_yticks(self.y_ticks)

        self.ax.hist(self.vel_mag_array[0], self.bins, edgecolor = 'black')
        
        self.anim = FuncAnimation(fig, self.update_graph, frames=self.vel_mag_array.shape[0],
                                  repeat=False, interval=10, blit=False)
        plt.show()
                
    
    def update_graph(self, i):
        '''
        used to update graph for FuncAnimation in animate_graph function.

        Parameters
        ----------
        i : int
            current time step

        Returns
        -------
        None
        '''
        self.ax.clear()
        self.ax.set_title('Velocity Distribution')
        self.ax.set_xlabel('Particle Velocity [m/s]')
        self.ax.set_ylabel('Number of Particles')
        self.ax.set_yticks(self.y_ticks)
        self.ax.hist(self.vel_mag_array[i], self.bins, edgecolor = 'black')
        self.ax.set_xlim(0, self.max_velocity + 10)
        self.ax.set_ylim(0, 120)
        
    def kinetic_energy_graph(self):
        '''
        plots total kinetic energy against time

        Returns
        -------
        None
        '''
        
        particle_mass = 1e-26
        ke_particles = 1/2 * particle_mass * (self.vel_mag_array**2)

        ''' sum all energies for per timestep'''
        ke_total = np.sum(ke_particles, axis = 1)
        
        plt.plot(ke_total)
        plt.xlabel('Timestep')
        plt.ylabel('Energy (J)')
                
        
        
        
        
        
        
if __name__ == "__main__":
    
    '''
    s = Simulation(1000, 0.02, 0.00002, 0.015, 7)
    s.simulate()
    
    '''
    data = np.load('MB_data_1000.npz')

    pos_array = data['pos']
    vel_array = data['vel']
    n = data['n']
    radius = data['r']
    box_size = data['box_size']
        
    
    animation = Animation(pos_array, vel_array, n, radius, box_size)
    animation.kinetic_energy_graph()
    

        
        
        
        
        
        
        

