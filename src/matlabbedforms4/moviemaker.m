% MOVIEMAKER - A subroutine of NEWDUNES that allows the user to determine the number and style of movie frames created
% File written by Carissa Carter 
% Last modified in 4/2005


NumberOfFrames = 1;
dT = 0;
dTrend = 0;
dZHO = 0;

switch colormode
    case {'color','c','Color','COLOR'}
% %************************************
% %MOVIE CONTROLS BY THE USER
k = input('How many frames showing deposition? :');    %number of depositional frames USER CHANGE
n = input('How many pause frames after deposition? :');    %length of pause after deposition  USER CHANGE
f = input('How many frames showing rotation? :');     %number of rotation frames USER CHANGE
r = input('How much total rotation in degrees? :');     %rotationin degrees  USER CHANGE
t = input('How many pause frames after rotation? :');    %length of pause after rotation  USER CHANGE
w = input('How many frames showing erosion? :');    %number of erosion frames  USER CHANGE
ii = input('How much erosion between frames? (0.025 suggested) :');    %increment of erosion USER CHANGE
d = input('How many pause frames at the end? :');
    case {'bw','BW','Bw','postscript','black and white'}
        k=1
        n=0
        f=0
        r=0
        t=0
        w=0
        ii=0
        d=0
end
%don't need to adjust this part
m = k + n;   
q = m + f;
j = r/f;       %rotation increment calculated
s = q + t;
v = s + w;     
b = -ii*w;      %end erosion value calculated
p = v + d;  %total number of frames
switch colormode
    case {'color','c','Color','COLOR'}
    Total_Frames_Per_Movie=p
end
dT(1:k) = 1:k;
dT(k+1:m) = k;
dT(m+1:p) = k;

dTrend(1:m) = 0;
dTrend(m+1:q) = j:j:r;
dTrend(q+1:s) = r;
dTrend(s+1:p) = r;

dZHO(1:s) = 0;
dZHO(s+1:v) = -ii:-ii:b;
dZHO(v+1:p) = b;


NumberOfFrames = length(dT);

%************************************