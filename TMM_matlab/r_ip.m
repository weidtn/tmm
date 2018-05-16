function y=r_ip(N0,N1,phi0,phi1)
%
% calculates the interface matrix for parallel polarization for
% the interface between medium 0 (incidence) and medium 1
%
% input parameters: 	N0	refr. index medium 1
%			N1
%			phi0	diffraction angle
%			phi1	angle of incidence
%
% returned result is a (2x2)-matrix

% fresnell coefficients according Azzam

rp=( N1*cos(phi0) - N0*cos(phi1) ) / ( N1*cos(phi0) + N0*cos(phi1) );
tp=2*N0*cos(phi0) / (N1*cos(phi0) + N0*cos(phi1));

% interface matrix

y = 1/tp*[ 1 rp; rp 1 ];

