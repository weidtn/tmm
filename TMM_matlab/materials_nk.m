function [lambda, n_material, k_material]=materials_nk(lambda_min,lambda_step,lambda_max,nk_file)

% interpolate values of a ellipsometer ascii file (lambda, n, k)
% usage e.g. 
% nk_path=['d:\Simulationen\Matlab\ellipsometry\nk_a-si.dat'];
% [lam,nmat,kmat]=materials_nk(500,1,1600,nk_path);


material=load(nk_file,'-ascii')';

%For testing purpose:
%plot(material(1,:),material(2,:),'b',material(1,:),material(3,:),'r');


lambda=lambda_min:lambda_step:lambda_max;
n_material=interp1(material(1,:),material(2,:),lambda);
k_material=interp1(material(1,:),material(3,:),lambda);