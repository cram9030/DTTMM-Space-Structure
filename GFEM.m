%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dx = GFEM(t,x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GFEM: Gaklerkin Finite Element Method					
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This funiction is called by an Matlab ODE solver to solve the beam
% bending problem proposed by a simple space structure. Resulting in a
% free-free boundary condition.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Author:	Nick Cramer, UCSC, Department of Computer Engineering			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Called by:	ODE45, ODE2, or ODE3											
% Calls:	MATLAB 5.2 std fcns
% Inputs:   t - time
%           x - beam states, velocities stacked on position vectors
% Outputs:  dx - actually dx/dt the rate of change of the states in time or
%                the temporal derivatives of them
% Globals: Due to the means that the odes work it is easier to define the
%          matrices as global variables than pass them another way
%           K - Stiffness matrix
%           C - Structural damping matrix
%           M - Mass matrix

global A B_dist

[row,col] = size(B_dist);

if t<1e-5
    dx = A*x+B_dist*rand(col,1);
else
    dx = A*x;
end