#include <dolfin.h>
#include "HyperElasticity.h"

using namespace dolfin;

class Rotation : public Expression
{
public:

  Rotation() : Expression(3) {}

  void eval(Array<double>& values, const Array<double>& x) const
  {
    const double scale = 0.5;

    const double x0 = -45;
    const double y0 = -205;

    double theta = M_PI * 10 / 180;

    double X = x0 + (x[0] - x0)*cos(theta) - (x[1] - y0)*sin(theta);
    double Y = y0 + (x[0] - x0)*sin(theta) + (x[1] - y0)*cos(theta);

    values[0] = scale*(X - x[0]);
    values[1] = scale*(Y - x[1]);
    values[2] = 0.0;
  }
};

class Bottom : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    return (x[2]  < 0) ;
  }
};
class Top : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    return (x[2] > 50) ;
  }
};

class Clamp : public Expression
{
public:

  Clamp() : Expression(3) {}

  void eval(Array<double>& values, const Array<double>& x) const
  {
    values[0] = 0.0;
    values[1] = 0.0;
    values[2] = 0.0;
  }

};
class BodyForce : public Expression
{
public:

  BodyForce() : Expression(3) {}

  void eval(Array<double>& values, const Array<double>& x) const
  {
    values[0] = 0.5 * cos(x[2]);
    values[1] = 0.5 * cos(x[2]);
    values[2] = sin(x[2]);
  }

};

int main()
{

  // Create mesh and define function space
  auto mesh = std::make_shared<Mesh>();
  auto mesh_file = std::make_shared<XDMFFile>(MPI_COMM_WORLD, "mesh/eevee.xdmf");
  mesh_file->read(*mesh);

  auto V = std::make_shared<HyperElasticity::FunctionSpace>(mesh);

  auto bottom = std::make_shared<Bottom>();
  auto top = std::make_shared<Top>();
  auto c = std::make_shared<Clamp>();
  auto r = std::make_shared<Rotation>();

  DirichletBC bcb(V, r, bottom);
  DirichletBC bct(V, c, top);
  std::vector<const DirichletBC*> bcs = {{&bcb, &bct}};


  auto B = std::make_shared<BodyForce>();
  auto T = std::make_shared<Constant>(0.0,  0.0, 0.0);

  auto u = std::make_shared<Function>(V);

  const double E  = 10.0;
  const double nu = 0.3;
  auto mu = std::make_shared<Constant>(E/(2*(1 + nu)));
  auto lambda = std::make_shared<Constant>(E*nu/((1 + nu)*(1 - 2*nu)));

  HyperElasticity::ResidualForm F(V);
  F.mu = mu; F.lmbda = lambda; F.u = u;
  F.B = B; F.T = T;

  HyperElasticity::JacobianForm J(V, V);
  J.mu = mu; J.lmbda = lambda; J.u = u;

  solve(F == 0, *u, bcs, J);

  File file("displacement.pvd");
  file << *u;

  return 0;
}
