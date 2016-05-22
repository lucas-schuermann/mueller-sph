#if __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <iostream>
#include <vector>
using namespace std;

#include <eigen3/Eigen/Dense>
using namespace Eigen;

// "Particle-based Viscoelastic Fluid Simluation"
// solver parameters
const static Vector2d G(0.f, -.02f * .25f); // external (gravitational) forces
const static float spacing = 2.f; // particle spacing/radius
const static float k = spacing / 1000.f; // far pressure weight
const static float k_near = k*10.f; // near pressure weight
const static float rest_density = 3.f;	// rest density
const static float r = spacing*1.25f; // kernel radius
const static float rsq = r*r; // radius^2 for optimization
const static float SIGMA = 0.2f; // visc parameters
const static float BETA = 0.2f;

// simulation parameters
const static float EPS = 1.0f; // boundary epsilon
const static float SPRING_CONST = 1./8.;
const static float MAX_VEL = 2.f; // velocity limit for stability
const static float MAX_VEL_SQ = MAX_VEL*MAX_VEL;

// neighbor data structure
// stores index of neighbor particles along with dist and squared dist to particle
struct Neighbor { 
	int j;
	float q, q2;
};

// particle data structure
// stores position, old position, velocity, and force for Verlet integration
// stores mass, rho, rho_near, pressure, pressure_near, sigma, and beta values for SPH
// stores list of neighboring particles for quick access in multiple simulation steps
struct Particle {
	Particle(float _x, float _y) : x(_x,_y), x0(_x,_y), v(0.f,0.f), f(0.f,0.f), rho(0.f), rho_near(0.f), p(0.f), p_near(0.f) {}
	Vector2d x, x0, v, f;
	float mass, rho, rho_near, p, p_near;
	vector<Neighbor> neighbors;
};

// solver data
const static int MAX_PARTICLES = 1000;
static vector<Particle> particles;

// rendering projection parameters
const static int WINDOW_WIDTH = 1280;
const static int WINDOW_HEIGHT = 800;
const static double VIEW_WIDTH = 50.0f;
const static double VIEW_HEIGHT = WINDOW_HEIGHT*VIEW_WIDTH/WINDOW_WIDTH;

void InitSPH(void)
{
	for(float y = EPS; y < VIEW_HEIGHT-EPS; y += r*0.5f)
		for(float x = EPS; x <= VIEW_WIDTH/2; x += r*0.5f)
			particles.push_back(Particle(x,y));
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
		glVertex2f(p.x(0), p.x(1));
	glEnd();

	glutSwapBuffers();
}

// Verlet integration
// spring boundary handling and sanity check
void Integrate(void)
{
	for(auto &p : particles)
	{
		// for simplicity, dt=1 assumed in integration
		p.x0 = p.x;
		p.x += p.v;
		p.x += p.f;

		// apply external forces
		p.f = G;
		p.v = p.x - p.x0;

		// if the velocity is greater than the max velocity, then cut it in half
		if(p.v.squaredNorm() > MAX_VEL_SQ)
			p.v = p.v * .5f;

		// enforce boundary condition
		if(p.x(0)-EPS < 0.0f)
			p.f(0) -= (p.x(0)-EPS) * SPRING_CONST;
		if(p.x(0)+EPS > VIEW_WIDTH) 
			p.f(0) -= (p.x(0)+EPS - VIEW_WIDTH) * SPRING_CONST;
		if(p.x(1)-EPS < 0.0f)
			p.f(1) -= (p.x(1)-EPS) * SPRING_CONST;
		if(p.x(1)+EPS > VIEW_HEIGHT)
			p.f(1) -= (p.x(1)+EPS - VIEW_HEIGHT) * SPRING_CONST;
	}
}

// LVSTODO: better description
// Calculate the density by basically making a weighted sum
// of the distances of neighboring particles within the radius of support (r)
void DensityStep(void)
{
	for(int i = 0; i < particles.size(); i++)
	{
		Particle &pi = particles[i];
		pi.rho = pi.rho_near = 0;
        
		// We will sum up the 'near' and 'far' densities.
		float d = 0.f, dn = 0.f;
        
		// Now look at every other particle
		pi.neighbors.clear();
		for(int j = i + 1; j < particles.size(); j++)
		{
			Particle &pj = particles[j];

			Vector2d rij = pj.x - pi.x;
			float rij_len2 = rij.squaredNorm();
			if(rij_len2 < rsq)
			{
				float rij_len = sqrt(rij_len2);
                
                // weighted distance
				float q = 1.f - rij_len / r;
				float q2 = q*q;
				float q3 = q2*q;
                
				d += q2;
				dn += q3;
                
				pj.rho += q2;
				pj.rho_near += q3;
                
				// add to neighbor list for faster access later
				Neighbor n;
				n.j = j;
				n.q = q; 
				n.q2 = q2;
				pi.neighbors.push_back(n);
			}
		}
        
		pi.rho += d;
		pi.rho_near += dn;
	}
}

void PressureAndViscosity(void)
{
	// pressure
	for(auto &pi : particles)
	{
		pi.p = k * (pi.rho - rest_density);
		pi.p_near = k_near * pi.rho_near;

		Vector2d dX = Vector2d(0,0);
		for(auto &n : pi.neighbors)
		{
			Particle &pj = particles[n.j];

			Vector2d rij = pj.x - pi.x;
			float dm = (pi.p + pi.p) * n.q + (pi.p_near + pj.p_near) * n.q2;
            
			Vector2d D = rij.normalized() * dm;
			dX += D;
			pj.f += D;
		}
		pi.f -= dX;
	}
    
	// viscosity
	for(auto &pi : particles)
	{
		for(auto &n : pi.neighbors)
		{
			Particle &pj = particles[n.j];   

			Vector2d rij = pj.x - pi.x;
			float l = (rij).norm();
			float q = l / r;
            
			Vector2d rijn = (rij / l);
			// Get the projection of the velocities onto the vector between them.
			float u = (pi.v - pj.v).dot(rijn);
			if(u > 0.f)
			{
				// Calculate the viscosity impulse between the two particles
				// based on the quadratic function of projected length.
				Vector2d I = (1.f - q) * (SIGMA * u + BETA * u*u) * rijn;
                
				// Apply the impulses on the two particles
				pi.v -= I * 0.5f;
				pj.v += I * 0.5f;
			}
            
		}
	}
}

void Update(void)
{ 
 	Integrate();
 	DensityStep();
    PressureAndViscosity();

	glutPostRedisplay();
}

void Keyboard(unsigned char c, __attribute__((unused)) int x, __attribute__((unused)) int y)
{   
	switch(c)
	{
	case ' ':
		if(particles.size() >= MAX_PARTICLES)
			std::cout << "maximum number of particles reached" << std::endl;
		else
			for(float y = VIEW_HEIGHT/1.5f-VIEW_HEIGHT/5.f; y < VIEW_HEIGHT/1.5f+VIEW_HEIGHT/5.f; y += r*0.5f)
				for(float x = VIEW_WIDTH/2.f-VIEW_HEIGHT/5.f; x <= VIEW_WIDTH/2.f+VIEW_HEIGHT/5.f; x += r*0.5f)
					particles.push_back(Particle(x,y));
		break;
	}
}

int main(int argc, char** argv)
{
	glutInitWindowSize(WINDOW_WIDTH,WINDOW_HEIGHT);
	glutInit(&argc, argv);
	glutCreateWindow("Basic SPH");
	glutDisplayFunc(Render);
	glutIdleFunc(Update);
	glutKeyboardFunc(Keyboard);

	InitGL();
	InitSPH();

	glutMainLoop();
	return 0;
}