from dolfin import *

def boundary(x, on_boundary):
    return on_boundary

def mesh_load():
    mesh = Mesh()
    XDMFFile("mesh/space.xdmf").read(mesh)
    return mesh

def function_space_generation(mesh):
    P1 = FiniteElement("CG", "triangle", 1)
    element =  MixedElement([P1, P1])
    G = FunctionSpace(mesh, P1)
    V = FunctionSpace(mesh, element)
    return G, V

def boundary_conditions(V):
    b_func_real = Constant(0.0);
    b_func_img = Constant(0.0);
    bc_real = DirichletBC(V.sub(0), b_func_real, boundary)
    bc_img = DirichletBC(V.sub(1), b_func_img, boundary)
    return [bc_real, bc_img]

def normalize_inplace(u):
    L = (u[0] * u[0] + u[1] * u[1]) * dx
    norm = sqrt(assemble(L))
    u.assign(u / norm)

def initial_distribution_choice(V):
    u_0 = Expression(("x[1] * exp(-1000 * (pow(x[0] - 0.5, 2) + pow(x[1], 2)))","0"),degree=2)
    u_n = interpolate(u_0, V)
    normalize_inplace(u_n)
    return u_n

T = 1
num_steps = 50
dt = T / num_steps
m = 1
omega = 10

G, V = function_space_generation(mesh_load())

v1, v2 = TestFunctions(V)
u = Function(V)
u_n = initial_distribution_choice(V)

zero = Expression("0", degree=2)
harmonic = Expression("m * pow(omega, 2) / 2 * (pow(x[0], 2) + pow(x[1], 2))", degree=1, m=m, omega=omega)
hydrogen = Expression("- 1 / (0.01 * pow(x[0] * x[0] + x[1] * x[1], 0.5) + 0.0001)", degree=2)
potential = Function(G)
potential.interpolate(harmonic)

F = (1/(2*m)*(inner(grad(u[0]), grad(v1)) + inner(grad(u[1]), grad(v2))) * dt
               + potential * (u[0] * v1 + u[1] * v2) * dt +
               ((u[1] - u_n[1])* v1 - (u[0] - u_n[0])* v2))*dx
bcs = boundary_conditions(V)

file_p = File("schrod/potential-harmonic.pvd")
file_p << potential
for i in range(num_steps):
    solve(F == 0, u, bcs)

    normalize_inplace(u)
    u_n.assign(u)

    file_psi = File("schrod/psi-harmonic-"+str(i)+".pvd")
    file_psi << u
