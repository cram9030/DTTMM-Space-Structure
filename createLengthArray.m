%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L = createLengthArray(beam_length,ActLocs,SegNum)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% createLengthArray: Create Length Array for GFEM			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is called to create an section length array for the
% Galerkin Finite Element Method with specified actuation locations for a
% given beam length and number of segements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Author:	Nick Cramer, UCSC, Department of Computer Engineering			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Called by:	CreateK_M_C											
% Calls:	MATLAB 5.2 std fcns
% Inputs:   beam_length - Length of beam
%           ActLocs -   Array of actuator locations starting from closest to
%                       origin and increasing as array indices increase (m)
%           SegNum -    Required number of segments for the half beam
% Outputs:  L - Array of Segment Lengths (m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes: assuming symettry of beam actuator locations
%       there is also some sort of error when the beam_length is to short
%       commented out code is to go to full beam instead of half

%Intialize variables
%num = SegNum/2;
num = SegNum;
L = zeros(SegNum,1);
unit = beam_length/SegNum; %calculate base unit
seg_Num_Act = (min(ActLocs)-mod(min(ActLocs),unit))/unit; %find how many base units can fit in the segment before the actuation locations

%Check different possible states, 1)to many segments, 2) just right number,
%3) to few segments
if length(ActLocs)+seg_Num_Act>num
    %When ther are to many segments remove non tip and actuator segments
    offset = length(ActLocs)+seg_Num_Act-num;
    L(1) = beam_length-max(ActLocs);
    L(2:length(diff(ActLocs))+1) = diff(ActLocs);
    L(length(diff(ActLocs))+2:num) = min(ActLocs)/(seg_Num_Act-offset)*ones(seg_Num_Act-offset,1);
    %L(num+1:end) = flipud(L(1:num));
elseif length(ActLocs)+seg_Num_Act == num
    L(1) = beam_length-max(ActLocs);
    L(2:length(diff(ActLocs))+1) = flipud(diff(ActLocs));
    L(length(diff(ActLocs))+2) = (beam_length-sum(L(1:length(diff(ActLocs))+1))-unit*(seg_Num_Act-1));
    L(length(diff(ActLocs))+3:num) = unit*ones(seg_Num_Act-1,1);
    %L(num+1:end) = flipud(L(1:num));
else
    %When there are too few put in all unit segments first
    remainder = num-(length(ActLocs)+seg_Num_Act);
    L(remainder+length(ActLocs)+2:num) = unit*ones(seg_Num_Act-1,1);
    L(remainder+length(ActLocs)+1) = ActLocs(1)-sum(unit*ones(seg_Num_Act-1,1));
    tempL = [beam_length-ActLocs(end);abs(diff(flip(ActLocs)))];
    %Then loop through the largest remaining segments halfing them until
    %the number of segments is met
    while remainder > 0
        [M,I] = max(tempL);
        tempL = [tempL(1:I-1);tempL(I)/2;tempL(I)/2;tempL(I+1:end)];
        remainder = remainder-1;
    end
    L(1:length(tempL)) = tempL;
    %L(num+1:end) = flipud(L(1:num));
end