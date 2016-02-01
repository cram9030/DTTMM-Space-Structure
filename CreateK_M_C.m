%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K,M,C] = CreateK_M_C(SegNum,beam_length,ActLocs,I,A,E,rho,etaE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CreateK_M_C: Stiffness, Mass, and Damping Matrix			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is called to create an section length array for the
% Galerkin Finite Element Method with specified actuation locations for a
% given beam length and number of segements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Author:	Nick Cramer, UCSC, Department of Computer Engineering			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Called by:    ------Fill in Later-----									
% Calls:	MATLAB 5.2 std fcns
% Inputs:   beam_length - Length of beam
%           ActLocs -   Array of actuator locations starting from closest to
%                       origin and increasing as array indices increase
%           SegNum -    Required number of segments for the half of the beam
%           I -         Beam area moment of inertia (m^4)
%           A -         Cross sectional area (m^2)
%           E -         Modulus of elasticity (Pa)
%           rho -       Density (kg/m^3)
%           etaE -      Material damping coefficent (Pa/s)
% Outputs:  K - Stiffness matrix
%           M - Mass matrix
%           C - Damping Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes: Only use for 2 segments or more and even Segment Number

L = createLengthArray(beam_length/2,ActLocs,SegNum);

K = zeros(2*(SegNum+1));
M = zeros(2*(SegNum+1));
C = zeros(2*(SegNum+1));
ksub1 = zeros(12);
msub1 = zeros(12);
csub1 = zeros(12);
ksub2 = zeros(12);
msub2 = zeros(12);
csub2 = zeros(12);

for i = 1:SegNum
    ksub1 = E*I/L(i)^3*[12 6 -12 6;...
        6 4 -6 2;...
        -12 -6 12 -6;...
        6 2 -6 4;...
        ];
    csub1 = etaE*I/L(i)^3*[12 6 -12 6;...
        6 4 -6 2;...
        -12 -6 12 -6;...
        6 2 -6 4;...
        ];
    msub1 = rho*A/(420*L(i))*[156 22 54 -13;...
        22 4 13 -3;...
        54 13 156 -22;...
        -13 -3 -22 4];

    K(2*(i-1)+1:2*(i-1)+4,2*(i-1)+1:2*(i-1)+4)=...
        K(2*(i-1)+1:2*(i-1)+4,2*(i-1)+1:2*(i-1)+4)+ksub1;
    C(2*(i-1)+1:2*(i-1)+4,2*(i-1)+1:2*(i-1)+4)=...
        C(2*(i-1)+1:2*(i-1)+4,2*(i-1)+1:2*(i-1)+4)+csub1;
    M(2*(i-1)+1:2*(i-1)+4,2*(i-1)+1:2*(i-1)+4) =...
        M(2*(i-1)+1:2*(i-1)+4,2*(i-1)+1:2*(i-1)+4)+msub1;
end