 function y=r_layer(N,d,phi,lambda)
%
% calculates a layer matrix according to Azzam
%
% input parameters: 	N  	complex refr. index
%			d  	layer thickness
%			phi 	diffraction angle (degrees)
%			lambda	wavelength [nm]
%

a=(2*pi*d*N/lambda)*cos(phi);

y = [exp(i*a) 0 ; 0 exp(-i*a)];
