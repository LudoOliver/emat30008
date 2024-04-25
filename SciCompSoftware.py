"""
This module contains the code for EMAT30008
Included are:
>ODESolver:
    >Solve_to()
        An ODE integrator, with options for both Euler and RungeKutta4 methods
>Shooting():
    Preforms numerical shooting for a given ODE
>Continuation:
    >NaturalContinuation()
        Natural parameter continuation of a given ODE
    >ShootingArcLength()
        Pseudo arc-length continuation of a given ODE
>FiniteDifferences()
    Solves a given ODE BVP using finite difference methods
>PDEMethods()
    Solves a given PDE using explicit/implicit Euler or RungeKutta4 methods

Examples and instructions:
> The convention for all ODEs is the form f(x,t,parameters*)
    Only continuation methods evaluate these parameters.
    Lambda functions should be used to fix parameter values outside of continuation.
> Boundary Conditions
    Three types of boundary condition are available and should be input as follows
    > Robin: du/dx|a = delta - gamma*u(a)
        BC = ("R",(delta,gamma))
    > Neumann: du/dx|a = delta
        BC = ("N",delta)  
    > Dirichlet: u(a) = alpha
        BC = ("D",alpha)
    The location of the boundary condition should be specified separately


"""
from Week16General import Shooting
import Week17Continuation as Continuation
import ODESolver
from FiniteDifferences import FiniteDifferences
from AMatrixBVector import MakeAMatrixBVector
import ExplicitEuler as PDEMethods



