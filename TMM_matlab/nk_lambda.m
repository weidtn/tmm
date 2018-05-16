function [lambda, n_material, k_material]=nk_lambda(lambda,nk_file)

% interpolate one value of a ellipsometer ascii file (lambda, n, k)
% usage e.g. 
% nk_path=['d:\Simulationen\Matlab\ellipsometry\nk_a-si.dat'];
% [lam,nmat,kmat]=nk_lambda(1500,nk_path);


material=load(nk_file,'-ascii')';

n_material=interp1(material(1,:),material(2,:),lambda);
k_material=interp1(material(1,:),material(3,:),lambda);