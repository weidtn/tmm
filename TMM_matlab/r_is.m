function y=r_ip(N0,N1,phi0,phi1)
%
% calculates the interface matrix for perpend. polarization at
% the interface between medium 0 (incidence) and medium 1
%
% input parameters: 	N1	refr. index medium 1
%			N0
%			phi1	diffraction angle
%			phi0	angle of incidence
%
% returned result is a (2x2)-matrix

% fresnell coefficients according Azzam

rs=( N1*cos(phi0) - N0*cos(phi1) ) / ( N1*cos(phi0) + N0*cos(phi1) );
ts=2*N0*cos(phi0) / (N0*cos(phi0) + N1*cos(phi1));

% interface matrix

y = 1/ts*[ 1 rs rs 1 ];
