%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xL,xR,dxL,dxR,ddxL,ddxR,ddTheta,dTheta,Theta,e] = DTFETMMHou_Loop_SS_2015_1_30(dt,E,etaE,I,rho,CA,L,xL_0,dxL_0,ddxL_0,xR_0,dxR_0,ddxR_0,fn_L,fn_R,c_w,I_h,Tf,t_imp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DTFETMMHou_Loop_SS_2015_1_30: Space Structure Discrete Time Transfer Method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function solves a simulation of a Galerkin Finite element
% discrete time transfer matrix method space structure representation using
% the Hubolt integration scheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Author:	Nick Cramer, UCSC, Department of Computer Engineering			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Called by:							
% Calls:	MATLAB 5.2 std fcns
% Inputs:   dt - Time step
%           E -   Modulus of elasticity (Pa)
%           etaE = Damping coefficient (Pa/s)
%           I - Area moment of inertia (m^4)
%           rho - Density (kg/m^3)
%           CA - Cross sectional area (m^2)
%           L - Length array (m)
%           xL_0 - Inital position array for the left
%           dxL_0 - Inital velocity array for the left
%           ddxL_0 - Inital acceleration array for the left
%           xR_0 - Inital position array for the right
%           dxR_0 - Inital velocity array for the right
%           ddxR_0 - Inital acceleration array for the right
%           iter - Current iteration to know what version of Houbolt
%           integration to use
%           fn_L - Array of forces of the left
%           fn_L - Array of forces of the right
%           c_w - Center rigid body width
%           I_h - Moment of inertia for center rotation piece
%           Tf - Final time for simulation
%           t_imp - Amount of time for the intial impulse lasts
% Outputs:  xL - Current state vector for the left
%           xR - Current state vector for the right
%           dxL - Position derivative for the left
%           dxR - Position derivative for the right
%           ddxL - Position double derivative for the left
%           ddxR - Position double derivative for the right
%           ddTheta - Rigid body rotation double derivative
%           dTheta -  Rigid body rotation derivative
%           Theta -  Rigid body rotation
%           e - Elapse time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Intialize arrays
t = cputime;
T = 0:dt:Tf;
xL = zeros(4,length(xL_0)+1,length(T));
dxL = zeros(2,length(xL_0)+1,length(T));
ddxL = zeros(2,length(xL_0)+1,length(T));
xL(:,:,1) = [zeros(4,1),[xL_0;zeros(2,length(xL_0))]];
dxL(:,:,1) = [zeros(2,1),dxL_0];
ddxL(:,:,1) = [zeros(2,1),ddxL_0];
xR = zeros(4,length(xR_0)+1,length(T));
dxR = zeros(2,length(xR_0)+1,length(T));
ddxR = zeros(2,length(xR_0)+1,length(T));
xR(:,:,1) = [zeros(4,1),[xR_0;zeros(2,length(xR_0))]];
dxR(:,:,1) = [zeros(2,1),dxR_0];
ddxR(:,:,1) = [zeros(2,1),ddxR_0];
ddTheta = zeros(1,length(T));

%Loop through simulation time
for i = 1:length(T)
    %Check to see if impulse force should be applied
    if T(i) < t_imp
        %Check to see the amount of time steps that have passed then
        %perform single step simulation
        if i == 1
            xL(:,:,i+1) = DTTMMHou2016_1_30(dt,E,etaE,I,rho,CA,L,xL(1:2,2:end,1),dxL(:,2:end,1),ddxL(:,2:end,1),xL(1:2,2:end,1),ddxL(:,2:end,1),xL(1:2,2:end,1),i,fn_L);
            xR(:,:,i+1) = DTTMMHou2016_1_30(dt,E,etaE,I,rho,CA,L,xR(1:2,2:end,1),dxR(:,2:end,1),ddxR(:,2:end,1),xR(1:2,2:end,1),ddxR(:,2:end,1),xL(1:2,2:end,1),i,fn_R);
        elseif i == 2
            xL(:,:,i+1) = DTTMMHou2016_1_30(dt,E,etaE,I,rho,CA,L,xL(1:2,2:end,i),dxL(:,2:end,i),ddxL(:,2:end,i),xL(1:2,2:end,i-1),ddxL(:,2:end,i-1),xL(1:2,2:end,i-1),i,fn_L);
            xR(:,:,i+1) = DTTMMHou2016_1_30(dt,E,etaE,I,rho,CA,L,xR(1:2,2:end,i),dxR(:,2:end,i),ddxR(:,2:end,i),xL(1:2,2:end,i-1),ddxR(:,2:end,i-1),xR(1:2,2:end,i-1),i,fn_R);
        else
            xL(:,:,i+1) = DTTMMHou2016_1_30(dt,E,etaE,I,rho,CA,L,xL(1:2,2:end,i),dxL(:,2:end,i),ddxL(:,2:end,i),xL(1:2,2:end,i-1),ddxL(:,2:end,i-1),xL(1:2,2:end,i-2),i,fn_L);
            xR(:,:,i+1) = DTTMMHou2016_1_30(dt,E,etaE,I,rho,CA,L,xR(1:2,2:end,i),dxR(:,2:end,i),ddxR(:,2:end,i),xR(1:2,2:end,i-1),ddxR(:,2:end,i-1),xR(1:2,2:end,i-2),i,fn_R);
        end
    else
        %Check to see the amount of time steps that have passed then
        %perform the simulation
        if i == 1
            xL(:,:,i+1) = DTTMMHou2016_1_30(dt,E,etaE,I,rho,CA,L,xL(1:2,2:end,1),dxL(:,2:end,1),ddxL(:,2:end,1),xL(1:2,2:end,1),ddxL(:,2:end,1),xL(1:2,2:end,1),i,zeros(2,20));
            xR(:,:,i+1) = DTTMMHou2016_1_30(dt,E,etaE,I,rho,CA,L,xR(1:2,2:end,1),dxR(:,2:end,1),ddxR(:,2:end,1),xR(1:2,2:end,1),ddxR(:,2:end,1),xR(1:2,2:end,1),i,zeros(2,20));
        elseif i == 2
            xL(:,:,i+1) = DTTMMHou2016_1_30(dt,E,etaE,I,rho,CA,L,xL(1:2,2:end,i),dxL(:,2:end,i),ddxL(:,2:end,i),xL(1:2,2:end,i-1),ddxL(:,2:end,i-1),xL(1:2,2:end,i-1),i,zeros(2,20));
            xR(:,:,i+1) = DTTMMHou2016_1_30(dt,E,etaE,I,rho,CA,L,xR(1:2,2:end,i),dxR(:,2:end,i),ddxR(:,2:end,i),xR(1:2,2:end,i-1),ddxR(:,2:end,i-1),xR(1:2,2:end,i-1),i,zeros(2,20));
        else
            xL(:,:,i+1) = DTTMMHou2016_1_30(dt,E,etaE,I,rho,CA,L,xL(1:2,2:end,i),dxL(:,2:end,i),ddxL(:,2:end,i),xL(1:2,2:end,i-1),ddxL(:,2:end,i-1),xL(1:2,2:end,i-2),i,zeros(2,20));
            xR(:,:,i+1) = DTTMMHou2016_1_30(dt,E,etaE,I,rho,CA,L,xR(1:2,2:end,i),dxR(:,2:end,i),ddxR(:,2:end,i),xR(1:2,2:end,i-1),ddxR(:,2:end,i-1),xR(1:2,2:end,i-2),i,zeros(2,20));
        end
    end
    %Check to see the iteration number for the integration step
    if i == 1
        [AnL,BnL,DnL,EnL] = Houbolt(dt,xL(1:2,2:end,1),dxL(:,2:end,1),ddxL(:,2:end,1),xL(1:2,2:end,1),ddxL(:,2:end,1),xL(1:2,2:end,1),i);
        [AnR,BnR,DnR,EnR] = Houbolt(dt,xR(1:2,2:end,1),dxR(:,2:end,1),ddxR(:,2:end,1),xR(1:2,2:end,1),ddxR(:,2:end,1),xR(1:2,2:end,1),i);
    elseif i == 2
        [AnL,BnL,DnL,EnL] = Houbolt(dt,xL(1:2,2:end,i),dxL(:,2:end,i),ddxL(:,2:end,i),xL(1:2,2:end,i-1),ddxL(:,2:end,i-1),xL(1:2,2:end,i),i);
        [AnR,BnR,DnR,EnR] = Houbolt(dt,xR(1:2,2:end,i),dxR(:,2:end,i),ddxR(:,2:end,i),xR(1:2,2:end,i-1),ddxR(:,2:end,i-1),xR(1:2,2:end,i),i);
    else
        [AnL,BnL,DnL,EnL] = Houbolt(dt,xL(1:2,2:end,i),dxL(:,2:end,i),ddxL(:,2:end,i),xL(1:2,2:end,i-1),ddxL(:,2:end,i-1),xL(1:2,2:end,i-2),i);
        [AnR,BnR,DnR,EnR] = Houbolt(dt,xR(1:2,2:end,i),dxR(:,2:end,i),ddxR(:,2:end,i),xR(1:2,2:end,i-1),ddxR(:,2:end,i-1),xR(1:2,2:end,i-2),i);
    end
    
    %Complete integraton step
    dxL(:,:,i+1) = DnL*xL(1:2,:,i+1)+[[0;0],EnL];
    ddxL(:,:,i+1) = AnL*xL(1:2,:,i+1)+[[0;0],BnL];
    dxR(:,:,i+1) = DnR*xR(1:2,:,i+1)+[[0;0],EnR];
    ddxR(:,:,i+1) = AnR*xR(1:2,:,i+1)+[[0;0],BnR];
    ddTheta(i) = c_w/(2*I_h)*xR(3,i+1)+xR(4,i+1)-c_w/(2*I_h)*xL(3,i+1)-xL(4,i+1);
end

%Solve for rigid body states
dTheta = cumtrapz(T,ddTheta);
Theta = cumtrapz(T,dTheta);
e = cputime-t;