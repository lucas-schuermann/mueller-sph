#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <GLUT/glut.h>
using namespace std;

//////////////////////////////////////////
// A simple two dimensional vector class
struct Vec2 {
    float x,y;
    Vec2() :x(0),y(0) { }
    Vec2(float a, float b) : x(a), y(b) { }
    Vec2 operator+(const Vec2& b) const { return Vec2(x+b.x, y+b.y); }
    Vec2 operator-(const Vec2& b) const { return Vec2(x-b.x, y-b.y); }
    Vec2 & operator=(const Vec2& b) { x=b.x; y=b.y; return *this; }
    Vec2 & operator+=(const Vec2& b) { return *this = *this + b; }
    Vec2 & operator-=(const Vec2& b) { return *this = *this - b; }
    
    float operator*(const Vec2& b) const { return x*b.x + y*b.y; }
    Vec2 operator*(float b) const { return Vec2(x * b, y * b); }
    Vec2 operator/(float b) const { return Vec2(x / b, y / b); }
    float len2() const { return *this * *this; }
    float len() const { return sqrt(len2()); }
    Vec2 normal() const { return *this / len(); }
};
Vec2 operator*(float b, const Vec2& a) { return Vec2(a.x * b, a.y * b); }

//////////////////////////////////////////
// A structure for holding two neighboring particles and their weighted distances
struct neighbor { int i, j; float q, q2; };

// The particle structure holding all of the relevant information.
struct particle {
	Vec2 pos, pos_old, vel, force;
	float mass, rho, rho_near, press, press_near, sigma, beta;
	vector<neighbor> neighbors;
};
//////////////////////////////////////////
int window_w=512, window_h=512;				// Initial Size of the Window
int N = 1000;										// Number of Particles in the simulation

float G = .02f * .25f;						// Gravitational Constant for our simulation

float spacing = 2.f;							// Spacing of particles
float k = spacing / 1000.0f;			   // Far pressure weight
float k_near = k*10;							// Near pressure weight
float rest_density = 3;						// Rest Density
float r=spacing*1.25f;						// Radius of Support
float rsq=r*r;									// ... squared for performance stuff

float SIM_W=50;								// The size of the world
float bottom = 0;								// The floor of the world

// Our collection of particles
vector<particle> particles;

// Mouse attractor
Vec2 attractor(999,999);
bool attracting = false;

//////////////////////////////////////////
// Some utility stuff.
// Between [0,1]
float rand01() { return (float)rand() * (1.f / RAND_MAX); }
// Between [a,b]
float randab(float a, float b) { return a + (b-a)*rand01(); }
//////////////////////////////////////////

void render()
{
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glLoadIdentity();
    
	// Draw Fluid Particles
	glPointSize(r*2);
	glBegin(GL_POINTS);
	for(int i=0; i < particles.size(); ++i)
	{
        glColor3f(.2f,.3f,.7f);
		glVertex2f(particles[i].pos.x, particles[i].pos.y);
	}
	glEnd();
    
	glutSwapBuffers();
}
bool asd = false;
//////////////////////////////////////////
void idle()
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
		particles[i].force = Vec2(0,-G);
        
		// Calculate the velocity for later.
		particles[i].vel = particles[i].pos - particles[i].pos_old;
        
		// If the velocity is really high, we're going to cheat and cap it.
		// This will not damp all motion. It's not physically-based at all. Just
		// a little bit of a hack.
		float max_vel = 2.f;
		float vel_mag = particles[i].vel.len2();
		// If the velocity is greater than the max velocity, then cut it in half.
		if(vel_mag > max_vel*max_vel)
			particles[i].vel = particles[i].vel * .5f;
        
		// If the particle is outside the bounds of the world, then
		// Make a little spring force to push it back in.
		if(particles[i].pos.x < -SIM_W) particles[i].force.x -= (particles[i].pos.x - -SIM_W) / 8;
		if(particles[i].pos.x >  SIM_W) particles[i].force.x -= (particles[i].pos.x - SIM_W) / 8;
		if(particles[i].pos.y < bottom) particles[i].force.y -= (particles[i].pos.y - bottom) / 8;
		if(particles[i].pos.y > SIM_W*2)particles[i].force.y -= (particles[i].pos.y - SIM_W*2) / 8;
        
		// Handle the mouse attractor.
		// It's a simple spring based attraction to where the mouse is.
		float attr_dist2 = (particles[i].pos - attractor).len2();
		const float attr_l = SIM_W/4;
		if( attracting )
			if( attr_dist2 < attr_l*attr_l )
				particles[i].force -= (particles[i].pos - attractor) / 256;
        
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
			Vec2 rij = particles[j].pos - particles[i].pos;
            
			// Along with the squared distance between
			float rij_len2 = rij.len2();
            
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
				neighbor n;
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
		Vec2 dX = Vec2();
        
		// For each of the neighbors
		unsigned long ncount = particles[i].neighbors.size();
		for(int ni=0; ni < ncount; ++ni)
		{
			neighbor n = particles[i].neighbors[ni];
			int j = n.j;
			float q(n.q);
			float q2(n.q2);
            
			// The vector from particle i to particle j
			Vec2 rij = particles[j].pos - particles[i].pos;
            
			// calculate the force from the pressures calculated above
			float dm = (particles[i].press + particles[j].press) * q +
            (particles[i].press_near + particles[j].press_near) * q2;
            
			// Get the direction of the force
			Vec2 D = rij.normal() * dm;
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
			neighbor n = particles[i].neighbors[ni];
            
			Vec2 rij = particles[n.j].pos - particles[i].pos;
			float l = (rij).len();
			float q = l / r;
            
			Vec2 rijn = (rij / l);
			// Get the projection of the velocities onto the vector between them.
			float u = (particles[n.i].vel - particles[n.j].vel) * rijn;
			if(u > 0)
			{
				// Calculate the viscosity impulse between the two particles
				// based on the quadratic function of projected length.
				Vec2 I = (1 - q) * (particles[n.j].sigma * u + particles[n.j].beta * u*u) * rijn;
                
				// Apply the impulses on the two particles
				particles[n.i].vel -= I * 0.5f;
				particles[n.j].vel += I * 0.5f;
			}
            
		}
	}
    
	// Draw the scene
	render();
}

//////////////////////////////////////////
void keyboard(unsigned char c, __attribute__((unused)) int x, __attribute__((unused))  int y)
{
    float radius = SIM_W/8;
    
	switch(c)
	{
            // Quit
		case 27:
		case 'q':
		case 'Q':
			exit(0);
			break;
            
            // If we press the space key, add some particles.
        case ' ':
            for(float y=SIM_W*2 - radius; y <= SIM_W*2+radius; y+=r*.5f){
                for(float x=-radius; x <= radius; x+=r*.5f)
                {
                    particle p;
                    p.pos = p.pos_old = Vec2(x , y) + Vec2(rand01(), rand01());
                    p.force = Vec2(0,0);
                    
                    if( (p.pos - Vec2( 0, SIM_W*2 ) ).len2() < radius*radius )
                        particles.push_back(p);
                }
			}
			break;
	}
}

void motion(int x, int y)
{
	// This simply updates the location of the mouse attractor.
	float relx = (float)(x - window_w/2) / window_w;
	float rely = -(float)(y - window_h) / window_h;
	Vec2 mouse = Vec2(relx*SIM_W*2, rely*SIM_W*2);
	attractor = mouse;
}

void mouse(__attribute__((unused)) int button, int state, __attribute__((unused))  int x, __attribute__((unused))  int y)
{
	if(state == GLUT_DOWN) attracting = true;
	else
	{
		attracting = false;
		attractor = Vec2(SIM_W * 99, SIM_W * 99);
	}
}

//////////////////////////////////////////
void init()
{
	// create a world with dimensions x:[-SIM_W,SIM_W] and y:[0,SIM_W*2]
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(-SIM_W,SIM_W,0,2*SIM_W);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
    
	glPointSize(5.f);
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    
	// Initialize particles
	// We will make a block of particles with a total width of 1/4 of the screen.
    float w = SIM_W/2;
    for(float y=bottom+1; y <= 10000; y+=r*.5f)
        for(float x=-w; x <= w; x+=r*.5f)
        {
            if(particles.size() > N) break;
            
            particle p;
            p.pos = Vec2(x, y);
            p.pos_old = p.pos + 0.001f * Vec2(rand01(), rand01());
            p.force = Vec2(0,0);
            p.sigma = .2f;
            p.beta = 0.2f;
            particles.push_back(p);
        }
}

//////////////////////////////////////////
int main(int argc, char **argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA|GLUT_DOUBLE);
	glutInitWindowSize(window_w, window_h);
	glutCreateWindow("SPH");
    
	glutDisplayFunc(render);
	glutKeyboardFunc(keyboard);
	glutIdleFunc(idle);
	glutMotionFunc(motion);
	glutMouseFunc(mouse);
    
	init();
	glutMainLoop();
}