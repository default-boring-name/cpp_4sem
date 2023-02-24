#include <dolfin.h>
#include "TentVel.h"
#include "Pressure.h"
#include "Velocity.h"

using namespace dolfin;

class NoslipDomain : public SubDomain
{
    bool inside (const Array<double>& x, bool on_boundary) const
    {
        return (on_boundary && (x[0] > DOLFIN_EPS || x[0] < 1 - DOLFIN_EPS));
    }
};
class InflowDomain : public SubDomain
{
    bool inside (const Array<double>& x, bool on_boundary) const
    {
        return  x[0] < DOLFIN_EPS;
    }
};

class Force : public Expression
{
public:
    double t;
    double x0;
    double y0;
    Force(double x0, double y0) : Expression(2), t(0), x0(x0), y0(y0){};
    void eval (Array<double>& values, const Array<double>& x) const
    { 
        double dx = x[0] - x0;
        double dy = x[1] - y0;
        values[0] =  10 * exp(-10*(dx * dx + dy * dy)) * cos(30 * t) * copysign(1, dx);
        values[1] =  10 * exp(-10*(dx * dx + dy * dy)) * cos(30 * t) * copysign(1, dy);
    }
};
class InflowPressure : public Expression
{
public:
    double t;
    InflowPressure() : t(0){};
    void eval (Array<double>& values, const Array<double>& x) const
    { 
        values[0] =  exp(-t) * (cos(3 * (3 * x[1] - 5 * t)) +  cos(3 *( 3 *  x[1] + 5 * t)));
    }
};

int main(){
   auto mesh = std::make_shared<Mesh>();
   auto mesh_file = std::make_shared<XDMFFile>(MPI_COMM_WORLD, "mesh/yung.xdmf");
   mesh_file->read(*mesh);

   auto V = std::make_shared<Velocity::FunctionSpace>(mesh);
   auto Q = std::make_shared<Pressure::FunctionSpace>(mesh);

   double dt = 0.01;

   auto p_in = std::make_shared<InflowPressure>();
   auto zero = std::make_shared<Constant>(0.0);
   auto zero_vector = std::make_shared<Constant>(0.0, 0.0);

   auto noslip_domain = std::make_shared<NoslipDomain>();
   auto inflow_domain = std::make_shared<InflowDomain>();
    
   DirichletBC noslip(V, zero_vector, noslip_domain);
   DirichletBC inflow(Q, p_in, inflow_domain);
   std::vector<DirichletBC*> bcu = {&noslip};
   std::vector<DirichletBC*> bcp = {&inflow};

   auto u0 = std::make_shared<Function>(V);
   auto u1 = std::make_shared<Function>(V);
   auto p1 = std::make_shared<Function>(Q);

   auto k = std::make_shared<Constant>(dt);
   auto f = std::make_shared<Force>(3, 0);
   
   TentVel::BilinearForm a1(V, V);
   TentVel::LinearForm L1(V);
   Pressure::BilinearForm a2(Q, Q);
   Pressure::LinearForm L2(Q);
   Velocity::BilinearForm a3(V, V);
   Velocity::LinearForm L3(V);

   a1.k = k; L1.k = k; L1.u0 = u0; L1.f = f;
   L2.k = k; L2.u1 = u1;
   L3.k = k; L3.u1 = u1; L3.p1 = p1;

   Matrix A1, A2, A3;
   assemble(A1, a1);
   assemble(A2, a2);
   assemble(A3, a3);

   Vector b1, b2, b3;

   double t = dt;
   for (unsigned i = 0; i < 100; ++i){
       f->t = t;
       f->y0 = 0.5 * sin(20 * t);
       p_in->t = t;

       assemble(b1, L1);
       for (unsigned i = 0; i < bcu.size(); i++)
           bcu[i]->apply(A1, b1);
       
       solve(A1, *u1->vector(), b1, "gmres", "default");

       assemble(b2, L2);
       for (unsigned i = 0; i < bcp.size(); i++)
       {
           bcp[i]->apply(A2, b2);
           bcp[i]->apply(*p1->vector());
       }
       
       solve(A2, *p1->vector(), b2, "bicgstab", "default");

       assemble(b3, L3);
       for (unsigned i = 0; i < bcu.size(); i++)
           bcu[i]->apply(A3, b3);
       
       solve(A3, *u1->vector(), b3, "gmres", "default");

       File ufile("results/dynamics-" + std::to_string(i) + ".pvd");
       ufile << *u1;
       *u0 = *u1;
       t += dt;
   }

    return 0;
}
