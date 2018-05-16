% Script file to calculate the reflectivity of an arbitraty
% stack of layers

clear;

% start parameters
exact = 100;                % precision, e.g. 1/100th nm
tuning_wav = 1550*exact;	% center wavelength
xmin = 1000*exact;			% min wavelength
xmax = 2600*exact;          % max wavelength
step = 1;                   % step size (wavelength)
phi(1) = 0;                 % incident angle
periodstop = 3;             % numbers of periods top DBR
periodsbot = 3;             % numbers of periods bottom DBR

% preallocate Refl-vector -------------

  Refl = zeros((xmax-xmin)/step,1);

% index of reflection

   index(1)=1;			% air
   index(2)=3.167;		% InP

% absorption coefficients in [cm-1] !!!

	alpha(1) = 0;
	alpha(2) = 0;

%for varth  = th_min:th_step:th_max;		% start thickness loop

% complex index of reflection

   d(2)=tuning_wav/(4*index(2));        % InP Bragg
   d(3)=tuning_wav/(4*index(1));    	% Air Bragg
   d(4)=tuning_wav;                 	% Air Cavity

for lambda  = xmin:step:xmax;           % start wavelength loop

% complex index of refraction

	for k = 1:1:2
		N(k) = lambda*alpha(k)*10^(-7)/(4*pi);
		N(k) = index(k) - i*N(k);
	end

% angles of diffraction ( >>> first element phi(1) defined above )

	for k = 2:1:2
		phi(k)=index(1)*sin(phi(1))/index(k);
		phi(k)=asin(phi(k));
	end


% ---------------------S I M U L A T I O N ------------------

% interface matrices parallel polarization

	k=1; l=2 ; Ip12 = r_ip( N(k),N(l),phi(k),phi(l) ) ; % air/InP
	k=2; l=1 ; Ip21 = r_ip( N(k),N(l),phi(k),phi(l) ) ; % InP/air

% interface matrices for perpendicular polarization

%	k=1; l=2 ; Is12 = r_is( N(k),N(l),phi(k),phi(l) ) ;
%   ... to be finished if needed

% layer matrices

	k=2 ; L2 = r_layer( N(2), d(k), phi(2), lambda );   % InP Bragg
 	k=3 ; L3 = r_layer( N(1), d(k), phi(1), lambda );   % Air Bragg
 	k=4 ; L4 = r_layer( N(1), d(k), phi(1), lambda );   % Air Cavity

% scattering matrix

  Sp = Ip12*L2*(Ip21*L3*Ip12*L2)^(periodstop-1)*Ip21*L4*Ip12*(L2*Ip21*L3*Ip12)^(periodsbot-1)*L2*Ip21*L3*Ip12;

% reflection and transmission coefficients

	Rp = Sp(2,1) / Sp(1,1) ;
%	Tp = 1 / Sp(1,1) ;
%	Rs = Ss(2,1) / Ss(1,1) ;
%	Ts = 1 / Ss(1,1) ;


%    Reflectivity for the wavelength lamda

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
