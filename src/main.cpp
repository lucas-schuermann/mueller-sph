#if __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif

#include <iostream>
#include <vector>
#include <math.h>
using namespace std;

#include <eigen3/Eigen/Dense>
using namespace Eigen;

// neighbor data structure
struct Neighbor { 
	int i, j; 
	float q, q2;
};

// particle data structure
struct Particle {
	Vector2d pos, pos_old, vel, force;
	float mass, rho, rho_near, press, press_near, sigma, beta;
	vector<Neighbor> neighbors;
};

// solver parameters
const static float G = .02f * .25f; // gravitational constant
const static float spacing = 2.f; // particle spacing/radius
const static float k = spacing / 1000.f; // far pressure weight
const static float k_near = k*10.f; // near pressure weight
const static float rest_density = 3.f;	// rest density
const static float r = spacing*1.25f; // kernel radius
const static float rsq = r*r;

// solver data
const static int MAX_PARTICLES = 1000;
static vector<Particle> particles;

// rendering projection parameters
const static int WINDOW_WIDTH = 1280;
const static int WINDOW_HEIGHT = 800;
const static double VIEW_WIDTH = 100.0f;
const static double VIEW_HEIGHT = WINDOW_HEIGHT*VIEW_WIDTH/WINDOW_WIDTH;

// LVSTODO: cleanup these...
float rand01() { return (float)rand() * (1.f / RAND_MAX); }
float randab(float a, float b) { return a + (b-a)*rand01(); }

void InitSPH(void)
{
	// Initialize particles
	// We will make a block of particles with a total width of 1/4 of the screen.
    float w = VIEW_WIDTH/2;
    for(float y=0.0f+1; y <= VIEW_HEIGHT; y+=r*.5f)
        for(float x=0.0; x <= w; x+=r*.5f)
        {
            if(particles.size() > MAX_PARTICLES) 
            	break;
            
            Particle p;
            p.pos = Vector2d(x, y);
            //p.pos_old = p.pos + 0.001f * Vector2d(rand01(), rand01());
            p.pos_old = p.pos;
            p.force = Vector2d(0,0);
            p.sigma = .2f;
            p.beta = 0.2f;
            particles.push_back(p);
        }
}

void InitGL(void)
{
	glClearColor(0.9f,0.9f,0.9f,1);
	glEnable(GL_POINT_SMOOTH);
	glPointSize(1.0f*WINDOW_WIDTH/VIEW_WIDTH);
	glMatrixMode(GL_PROJECTION);
}

void Render(void)
{
	glClear(GL_COLOR_BUFFER_BIT);
	
	glLoadIdentity();
	glOrtho(0, VIEW_WIDTH, 0, VIEW_HEIGHT, 0, 1);


	glColor4f(0.2f, 0.6f, 1.0f, 1);
	glBegin(GL_POINTS);
	for(auto &p : particles)
		glVertex2f(p.pos(0), p.pos(1));
	glEnd();

	glutSwapBuffers();
}

void Update(void)
{
	// UPDATE
	//
	// This modified verlet integrator has dt = 1 and calculates the velocity
	// For later use in the simulation.
    
	// For each particles i ...
    for(int i=0; i < particles.size(); ++i)
	{
		// Normal verlet stuff
		particles[i].pos_old = particles[i].pos;
		particles[i].pos += particles[i].vel;
        
		// Apply the currently accumulated forces
		particles[i].pos += particles[i].force;
        
		// Restart the forces with gravity only. We'll add the rest later.
		particles[i].force = Vector2d(0,-G);
        
		// Calculate the velocity for later.
		particles[i].vel = particles[i].pos - particles[i].pos_old;
        
		// If the velocity is really high, we're going to cheat and cap it.
		// This will not damp all motion. It's not physically-based at all. Just
		// a little bit of a hack.
		float max_vel = 2.f;
		float vel_mag = particles[i].vel.squaredNorm();
		// If the velocity is greater than the max velocity, then cut it in half.
		if(vel_mag > max_vel*max_vel)
			particles[i].vel = particles[i].vel * .5f;
        
		// If the particle is outside the bounds of the world, then
		// Make a little spring force to push it back in.
		float eps = 1.0f;
		if(particles[i].pos(0)-eps < 0.0f) particles[i].force(0) -= (particles[i].pos(0)-eps - 0.0f) / 8;
		if(particles[i].pos(0)+eps >  VIEW_WIDTH) particles[i].force(0) -= (particles[i].pos(0)+eps - VIEW_WIDTH) / 8;
		if(particles[i].pos(1)-eps < 0.0f) particles[i].force(1) -= (particles[i].pos(1)-eps - 0.0f) / 8;
		if(particles[i].pos(1)+eps > VIEW_HEIGHT)particles[i].force(1) -= (particles[i].pos(1)+eps - VIEW_HEIGHT) / 8;
        
        
		// Reset the nessecary items.
		particles[i].rho = particles[i].rho_near = 0;
		particles[i].neighbors.clear();
	}
    
	// DENSITY
	//
	// Calculate the density by basically making a weighted sum
	// of the distances of neighboring particles within the radius of support (r)
    
	// For each particle ...
	for(int i=0; i < particles.size(); ++i)
	{
		particles[i].rho = particles[i].rho_near = 0;
        
		// We will sum up the 'near' and 'far' densities.
		float d=0, dn=0;
        
		// Now look at every other particle
        unsigned long particles_size = particles.size();
		for(int j = i + 1; j < particles_size; ++j)
		{
			// The vector seperating the two particles
			Vector2d rij = particles[j].pos - particles[i].pos;
            
			// Along with the squared distance between
			float rij_len2 = rij.squaredNorm();
            
			// If they're within the radius of support ...
			if(rij_len2 < rsq)
			{
				// Get the actual distance from the squared distance.
				float rij_len = sqrt(rij_len2);
                
				// And calculated the weighted distance values
				float q = 1 - rij_len / r;
				float q2 = q*q;
				float q3 = q2*q;
                
				d += q2;
				dn += q3;
                
				// Accumulate on the neighbor
				particles[j].rho += q2;
				particles[j].rho_near += q3;
                
				// Set up the neighbor list for faster access later.
				Neighbor n;
				n.i = i; n.j = j;
				n.q = q; n.q2 = q2;
				particles[i].neighbors.push_back(n);
			}
		}
        
		particles[i].rho        += d;
		particles[i].rho_near   += dn;
	}
    
	// PRESSURE
	//
	// Make the simple pressure calculation from the equation of state.
	for(int i=0; i < particles.size(); ++i)
	{
		particles[i].press = k * (particles[i].rho - rest_density);
		particles[i].press_near = k_near * particles[i].rho_near;
	}
    
	// PRESSURE FORCE
	//
	// We will force particles in or out from their neighbors
	// based on their difference from the rest density.
    
	// For each particle ...
	for(int i=0; i < particles.size(); ++i)
	{
		Vector2d dX = Vector2d(0,0);
        
		// For each of the neighbors
		unsigned long ncount = particles[i].neighbors.size();
		for(int ni=0; ni < ncount; ++ni)
		{
			Neighbor n = particles[i].neighbors[ni];
			int j = n.j;
			float q(n.q);
			float q2(n.q2);
            
			// The vector from particle i to particle j
			Vector2d rij = particles[j].pos - particles[i].pos;
            
			// calculate the force from the pressures calculated above
			float dm = (particles[i].press + particles[j].press) * q +
            (particles[i].press_near + particles[j].press_near) * q2;
            
			// Get the direction of the force
			Vector2d D = rij.normalized() * dm;
			dX += D;
			particles[j].force += D;
		}
        
		particles[i].force -= dX;
	}
    
	// VISCOSITY
	//
	// This simulation actually may look okay if you don't compute
	// the viscosity section. The effects of numerical damping and
	// surface tension will give a smooth appearance on their own.
	// Try it.
    
	// For each particle
	for(int i=0; i < particles.size(); ++i)
	{
		// For each of that particles neighbors
		for(int ni=0; ni < particles[i].neighbors.size(); ++ni)
		{
			Neighbor n = particles[i].neighbors[ni];
            
			Vector2d rij = particles[n.j].pos - particles[i].pos;
			float l = (rij).norm();
			float q = l / r;
            
			Vector2d rijn = (rij / l);
			// Get the projection of the velocities onto the vector between them.
			float u = (particles[n.i].vel - particles[n.j].vel).dot(rijn);
			if(u > 0)
			{
				// Calculate the viscosity impulse between the two particles
				// based on the quadratic function of projected length.
				Vector2d I = (1 - q) * (particles[n.j].sigma * u + particles[n.j].beta * u*u) * rijn;
                
				// Apply the impulses on the two particles
				particles[n.i].vel -= I * 0.5f;
				particles[n.j].vel += I * 0.5f;
			}
            
		}
	}

	glutPostRedisplay();
}

void Keyboard(unsigned char c, __attribute__((unused)) int x, __attribute__((unused)) int y)
{
	float radius = VIEW_WIDTH/8;
    
	switch(c)
	{
	case ' ':
        float w = VIEW_WIDTH/4;
	    for(float y=VIEW_HEIGHT/2+1; y <= VIEW_HEIGHT; y+=r*.5f)
	        for(float x=0.0+1.0; x <= w; x+=r*.5f)
	        {
	            Particle p;
	            p.pos = Vector2d(x, y);
	            //p.pos_old = p.pos + 0.001f * Vector2d(rand01(), rand01());
	            p.pos_old = p.pos;
	            p.force = Vector2d(0,0);
	            p.sigma = .2f;
	            p.beta = 0.2f;
	            particles.push_back(p);
	        }
		break;
	}
}

int main(int argc, char** argv)
{
	glutInitWindowSize(WINDOW_WIDTH,WINDOW_HEIGHT);
	glutInit(&argc, argv);
	glutCreateWindow("SPH");
	glutDisplayFunc(Render);
	glutIdleFunc(Update);
	glutKeyboardFunc(Keyboard);

	InitGL();
	InitSPH();

	glutMainLoop();
	return 0;
}