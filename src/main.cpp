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

// "Particle-Based Fluid Simulation for Interactive Applications"
// solver parameters
const static Vector2d G(0.f, -10.f); // external (gravitational) forces
const static float REST_DENS = 1.f; // rest density
const static float GAS_CONST = 1.f; // const for equation of state
const static float H = 1.5f; // kernel radius
const static float HSQ = H*H; // radius^2 for optimization
const static float MASS = 0.02f; // assume all particles have the same mass
const static float VISC = 6.5f; // viscosity constant
const static float DT = 0.001f; // integration timestep

// smoothing kernels defined in Müller and their gradients
const static float POLY6 = 315.f/(65.f*M_PI*pow(H, 9.f));
const static float SPIKY_GRAD = -45.f/(M_PI*pow(H, 6.f));
const static float VISC_LAP = 45.f/(M_PI*pow(H, 6.f));

// simulation parameters
const static float EPS = 1.0f; // boundary epsilon
const static float BOUND_DAMPING = -0.5f;

// particle data structure
// stores position, velocity, and force for integration
// stores density (rho) and presure values for SPH
struct Particle {
	Particle(float _x, float _y) : x(_x,_y), v(0.f,0.f), f(0.f,0.f), rho(0.f), p(0.f) {}
	Vector2d x, v, f;
	float rho, p;
};

// solver data
const static int MAX_PARTICLES = 5;
static vector<Particle> particles;

// rendering projection parameters
const static int WINDOW_WIDTH = 1280;
const static int WINDOW_HEIGHT = 800;
const static double VIEW_WIDTH = 50.0f;
const static double VIEW_HEIGHT = WINDOW_HEIGHT*VIEW_WIDTH/WINDOW_WIDTH;

void InitSPH(void)
{
	for(float y = EPS+2.f; y < VIEW_HEIGHT-EPS; y += H*0.5f)
		for(float x = EPS+2.f; x <= VIEW_WIDTH/2; x += H*0.5f)
			if(particles.size() < MAX_PARTICLES)
				particles.push_back(Particle(x,y));
}

void Integrate(void)
{
	for(auto &p : particles)
	{
		// forward euler
		p.v += DT*p.f/p.rho;
		p.x += DT*p.v;

		// enforce boundary conditions
		if(p.x(0)-EPS < 0.0f)
		{
			p.v(0) *= BOUND_DAMPING;
			p.x(0) = EPS;
		}
		if(p.x(0)+EPS > VIEW_WIDTH) 
		{
			p.v(0) *= BOUND_DAMPING;
			p.x(0) = VIEW_WIDTH-EPS;
		}
		if(p.x(1)-EPS < 0.0f)
		{
			p.v(1) *= BOUND_DAMPING;
			p.x(1) = EPS;
		}
		if(p.x(1)+EPS > VIEW_HEIGHT)
		{
			p.v(1) *= BOUND_DAMPING;
			p.x(1) = VIEW_HEIGHT-EPS;
		}
	}
}

void ComputeDensityPressure(void)
{
	for(auto &pi : particles)
	{
		pi.rho = 0.f;
		for(auto &pj : particles)
		{
			if(&pi == &pj) break;

			float r2 = (pj.x - pi.x).squaredNorm();
			if(r2 < HSQ)
			{
				// this computation is symmetric
				pi.rho += MASS*POLY6*pow(HSQ-r2, 3.f);
			}
		}
		//pi.p = GAS_CONST*(pi.rho - REST_DENS);
		pi.p = GAS_CONST*(pow(pi.rho/REST_DENS,3.f)-1.f);
		//std::cout << "DENSITY" << std::endl;
		//std::cout << pi.rho << std::endl;
		//std::cout << pi.p << std::endl;
	}
}

void ComputeForces(void)
{
	// pressure
	for(auto &pi : particles)
	{
		Vector2d fpress(0.f, 0.f);
		Vector2d fvisc(0.f, 0.f);
		Vector2d fgrav(0.f, 0.f);
		for(auto &pj : particles)
		{
			if(&pi == &pj) break;

			float r = (pj.x - pi.x).norm();
			if(r < H)
			{
				// pressure
				fpress += -(pj.x - pi.x).normalized()*MASS*(pi.p + pj.p)/(2.f * pj.rho) * SPIKY_GRAD*pow(H-r,2.f);
				// visc
				fvisc += VISC*MASS*(pj.v - pi.v)/pj.rho * VISC_LAP*(H-r);
			}
		}
		std::cout << "FORCES" << std::endl;
		std::cout << fpress << std::endl;
		std::cout << fvisc << std::endl;
		fgrav = G * pi.rho;
		pi.f = fpress + fvisc + fgrav;
		std::cout << pi.f << std::endl;
	}
}

static bool once = false;
void Update(void)
{ 
	if(!once)
	{
 		ComputeDensityPressure();
   	 	ComputeForces();
   	 	Integrate();
   	 	once = false;
	}

	glutPostRedisplay();
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

void Keyboard(unsigned char c, __attribute__((unused)) int x, __attribute__((unused)) int y)
{   
	switch(c)
	{
	case ' ':
		if(particles.size() >= MAX_PARTICLES)
			std::cout << "maximum number of particles reached" << std::endl;
		else
			for(float y = VIEW_HEIGHT/1.5f-VIEW_HEIGHT/5.f; y < VIEW_HEIGHT/1.5f+VIEW_HEIGHT/5.f; y += H*0.5f)
				for(float x = VIEW_WIDTH/2.f-VIEW_HEIGHT/5.f; x <= VIEW_WIDTH/2.f+VIEW_HEIGHT/5.f; x += H*0.5f)
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