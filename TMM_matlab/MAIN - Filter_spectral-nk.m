% Script file to calculate the reflectivity of an arbitraty
% stack of layers

clear;

% start parameters

exact = 100;                % precision, e.g. 1/100th nm
tuning_wav = 1550*exact;	% center wavelength
xmin = 1000*exact;			% min wavelength
xmax = 1700*exact;          % max wavelength
step = 0.05*exact;          % step size (wavelength, e.g. in nm)
phi(1) = 0;                 % incident angle
periodstop = 4;             % numbers of periods top DBR
periodsbot = 4;             % numbers of periods bottom DBR
nk_path1=['d:\Simulationen\Matlab\ellipsometry\nk_sio2.dat'];  % path of ellipsometry file SiO2
nk_path2=['d:\Simulationen\Matlab\ellipsometry\nk_a-si.dat'];  % path of ellipsometry file a-Si

% preallocate Refl-vector -------------

  Refl = zeros((xmax-xmin)/step,1);

% interpolate n and k for simulated wavelengths

  [lam1,nmat1,kmat1]=materials_nk(xmin/exact,step/exact,xmax/exact,nk_path1);     % material L
  [lam2,nmat2,kmat2]=materials_nk(xmin/exact,step/exact,xmax/exact,nk_path2);     % material H

% interpolate n and k for reference wavelength

  [ref_lam1,ref_nmat1,ref_kmat1]=nk_lambda(tuning_wav/exact,nk_path1); % reference L
  [ref_lam2,ref_nmat2,ref_kmat2]=nk_lambda(tuning_wav/exact,nk_path2); % reference H

%for varth  = th_min:th_step:th_max;		% start thickness loop

% layer thicknesses

   d(2)=tuning_wav/(4*ref_nmat2);        % H QWOT
   d(3)=tuning_wav/(4*ref_nmat1);    	 % L QWOT
   d(4)=tuning_wav/(2*ref_nmat1);        % L Cavity

lix=0;      % loop counter

for lambda  = xmin:step:xmax;            % start wavelength loop

% complex index of refraction

    lix=lix+1;

	N(1) = nmat1(lix) - i*kmat1(lix);
    N(2) = nmat2(lix) - i*kmat2(lix);

% angles of refraction ( >>> first element phi(1) defined above )

    phi(2)=asin(nmat1(lix)*sin(phi(1))/nmat2(lix));


% ---------------------S I M U L A T I O N ------------------

% interface matrices for parallel polarization

	k=1; l=2 ; Ip12 = r_ip( N(k),N(l),phi(k),phi(l) ) ; % L/H
	k=2; l=1 ; Ip21 = r_ip( N(k),N(l),phi(k),phi(l) ) ; % H/L

% interface matrices for perpendicular polarization

%	k=1; l=2 ; Is12 = r_is( N(k),N(l),phi(k),phi(l) ) ;
%   ... to be finished if needed

% layer matrices

	k=2 ; L2 = r_layer( N(2), d(k), phi(2), lambda );   % H QWOT
 	k=3 ; L3 = r_layer( N(1), d(k), phi(1), lambda );   % L QWOT
 	k=4 ; L4 = r_layer( N(1), d(k), phi(1), lambda );   % L Cavity

% scattering matrix

  Sp = Ip12*L2*(Ip21*L3*Ip12*L2)^(periodstop-1)*Ip21*L4*Ip12*(L2*Ip21*L3*Ip12)^(periodsbot-1)*L2*Ip21*L3*Ip12;

% reflection and transmission coefficients

	Rp = Sp(2,1) / Sp(1,1) ;
%	Tp = 1 / Sp(1,1) ;
%	Rs = Ss(2,1) / Ss(1,1) ;
%	Ts = 1 / Ss(1,1) ;


%    Reflectivity for the wavelength lambda

	pos = (lambda-xmin+step)/step ;
    Refl(pos,1) = lambda/exact ;
    Refl(pos,2) = 100*(1-abs((Rp)^2)) ;


end					% end wavelength loop


%    Plot the result

   plot(Refl(:,1),Refl(:,2),'-')

%end

   xmin=xmin/exact;
   xmax=xmax/exact;
   xlabel('Wavelength [nm]')
   ylabel('Transmission [%]')
   axis([xmin xmax 0 105])
   title('Air Gap Filter')
