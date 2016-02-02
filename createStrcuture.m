%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,B_dist,B_cont,C,K_R,K_L,M_R,M_L] = createStrcuture(SegNum,beam_length,Locs,ControlLocs,DisturbLocs,SensorLocs,I,A,E,rho,etaE,c_w,I_h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% createStrcuture: Create GFEM space structure with rigid section			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function creates the the state, input, and output matrices of an space
% structure with a rigid center component. The mass and stiffness matrices
% are also returned for both the right and the left side of the space
% structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Author:	Nick Cramer, UCSC, Department of Computer Engineering			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Called by:    --------------------										
% Calls:	MATLAB 5.2 std fcns
%           CreateK_M_C
%           createLengthArray
% Inputs:   beam_length - Length of beam
%           Locs -   Array of actuator and sensor locations starting from closest to
%                       origin and increasing as array indices increase (m)
%           SegNum -    Required number of segments for the half beam
%           ControlLocs - Array of control input locations
%           DisturbLocs - Array of disturbance input locations
%           SensorLocs - Array of sensor locations
%           I - Area moment of inertia of flexible beam
%           A - Cross section area of flexible beam
%           E - Modulus of elasticity
%           rho - Density of beams
%           etaE - Damping coefficeint of beams
%           c_w - Center rigid component width
%           I_h - Moment of inertia for center
% Outputs:  A - State matrix of size 2*(4*SegNum+1)
%           B_dist - Disturbance input matrix of size 2*(4*SegNum+1) by
%                   length(DisturbLocs)
%           B_cont - Control input matrix of size 2*(4*SegNum+1) by
%                   length(ControlLocs)
%           C - Output matrix
%           K_R - Right stiffness matrix
%           K_L - Left stiffness matrix
%           M_R - Right mass matrix
%           M_L - Left mass matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Notes: Currently assumes symetric structure

%Intialize Matrices
AsubK = zeros(4*SegNum+1);
AsubC = zeros(4*SegNum+1);
controlInput = zeros(2*SegNum+2,length(ControlLocs));
disturbInupt = zeros(2*SegNum+2,length(DisturbLocs));
sensorOutput = zeros(2*SegNum+2,length(SensorLocs));

%Create left and right beam
[K_L,M_L,C_L] = CreateK_M_C(SegNum,beam_length,Locs,I,A,E,rho,etaE);
tempK_R = rot90(K_L,2);
K_R = tempK_R;
K_R(:,1:2:end) = tempK_R(:,2:2:end);
K_R(:,2:2:end) = tempK_R(:,1:2:end);
tempK_R = K_R;
K_R(1:2:end,:) = tempK_R(2:2:end,:);
K_R(2:2:end,:) = tempK_R(1:2:end,:);
tempM_R = rot90(M_L,2);
M_R = tempM_R;
M_R(:,1:2:end) = tempM_R(:,2:2:end);
M_R(:,2:2:end) = tempM_R(:,1:2:end);
tempM_R = M_R;
M_R(1:2:end,:) = tempM_R(2:2:end,:);
M_R(2:2:end,:) = tempM_R(1:2:end,:);
tempC_R = rot90(C_L,2);
C_R = tempC_R;
C_R(:,1:2:end) = tempC_R(:,2:2:end);
C_R(:,2:2:end) = tempC_R(:,1:2:end);
tempC_R = C_R;
C_R(1:2:end,:) = tempC_R(2:2:end,:);
C_R(2:2:end,:) = tempC_R(1:2:end,:);

%Enforce boundary conditions and add rigid body mechanics to the stiffness
%matrices
AsubK_L = -M_L\K_L;
AsubK_L = AsubK_L(1:2*SegNum,1:2*SegNum);
AsubK_R = -M_R\K_R;
AsubK_R = AsubK_R(3:end,3:end);
AsubK(1:2*SegNum,1:2*SegNum) = AsubK_L;
AsubK(end-2*SegNum+1:end,end-2*SegNum+1:end) = AsubK_R;
AsubK(2*SegNum+1,2*SegNum-1:2*SegNum) = -c_w/(2*I_h)*K_L(2*SegNum+1,2*SegNum-1:2*SegNum)-K_L(2*SegNum+2,2*SegNum-1:2*SegNum);
AsubK(2*SegNum+1,2*SegNum+2:2*SegNum+3) = c_w/(2*I_h)*K_R(1,3:4)+K_R(2,3:4);

%Enforce boundary conditions and add rigid body mechanics to the damping
%matrices
AsubC_L = -M_L\C_L;
AsubC_L = AsubC_L(1:2*SegNum,1:2*SegNum);
AsubC_R = -M_R\C_R;
AsubC_R = AsubC_R(3:end,3:end);
AsubC(1:2*SegNum,1:2*SegNum) = AsubC_L;
AsubC(end-2*SegNum+1:end,end-2*SegNum+1:end) = AsubC_R;
AsubC(2*SegNum+1,2*SegNum-1:2*SegNum) = -c_w/(2*I_h)*C_L(2*SegNum+1,2*SegNum-1:2*SegNum)-C_L(2*SegNum+2,2*SegNum-1:2*SegNum);
AsubC(2*SegNum+1,2*SegNum+2:2*SegNum+3) = c_w/(2*I_h)*C_R(1,3:4)+C_R(2,3:4);

%Create A matrix by combinging sub matrices
A = [zeros(size(AsubK)),eye(size(AsubK));AsubK,AsubC];

%Create input matrices for the disturbance and control locations prior to
%state stacking
L = createLengthArray(beam_length/2,Locs,SegNum);

for i = 1:length(DisturbLocs)
    segment = find(flipud(cumsum(flipud(L)))==DisturbLocs(i));
    disturbInupt(2*segment+2,i) = 1;
end
for i = 1:length(ControlLocs)
    segment = find(flipud(cumsum(flipud(L)))==ControlLocs(i));
    controlInput(2*segment-3,i) = 1;
    controlInput(2*segment-1,i) = -1;
end

%Invert mass matrix and stack states to create controlable input and
%disturbance matrices
%-------------------------------------------------------------------
disturbInupt = -M_L\disturbInupt;
%Create opposite side disturbance input
disturbInput_R = flipud(disturbInupt);
tempDisturbInput_R = disturbInput_R;
disturbInput_R(2:2:end,:) = tempDisturbInput_R(1:2:end,:);
disturbInput_R(1:2:end,:) = tempDisturbInput_R(2:2:end,:);

controlInput = -M_L\controlInput;
%Create opposite side control input
controlInput_R = flipud(controlInput);
tempControlInput_R = controlInput_R;
controlInput_R(2:2:end,:) = tempControlInput_R(1:2:end,:);
controlInput_R(1:2:end,:) = tempControlInput_R(2:2:end,:);

%Create B_disturbance and B_control Matrices
B_dist = [[disturbInupt(1:end-2,:),zeros(size(disturbInupt(1:end-2,:)))];zeros(1,2*length(DisturbLocs));[zeros(size(disturbInupt(1:end-2,:))),disturbInput_R(3:end,:)]];
B_dist = [zeros(size(B_dist));B_dist];
B_cont = [[controlInput(1:end-2,:),zeros(size(controlInput(1:end-2,:)))];zeros(1,2*length(ControlLocs));[zeros(size(controlInput(1:end-2,:))),controlInput_R(3:end,:)]];
B_cont = [zeros(size(B_cont));B_cont];

%Create Output Matrix assuming that the sensors being used are only rate
%gyros
for i = 1:length(SensorLocs)
    segment = find(flipud(cumsum(flipud(L)))==SensorLocs(i));
    sensorOutput(2*segment,i) = 1;
end
%Create opposite side output sub matrix
sensorOutput_R = flipud(sensorOutput);
tempSensorOutput_R = sensorOutput_R;
sensorOutput_R(2:2:end,:) = tempSensorOutput_R(1:2:end,:);
sensorOutput_R(1:2:end,:) = tempSensorOutput_R(2:2:end,:);

%Assemble C output matrix
[row,col] = size(sensorOutput(1:end-2,:));
C = [[sensorOutput(1:end-2,:),zeros(row,col+2)];[zeros(1,length(SensorLocs)),1,0,zeros(1,length(SensorLocs))];[zeros(row,col+2),sensorOutput_R(3:end,:)]]';
C = [zeros(size(C)),C];
C(length(SensorLocs)+1,2*SegNum+1) = 1;