tic

% This file will calculate the measurement signal and minimum resolvable
% displacement for a cubic particle in optical tweezers, with particle
% tracking via a quadrant detector at the back-focal plane of a condenser.

% This file is part of the Quadrant Detection Package version 1.1. 
% M. A. Taylor and W. P. Bowen, “A computational tool to characterize
% particle tracking measurements in optical tweezers,” J. Opt. 15, 085701
% (2013).  
% https://github.com/michael-a-taylor/QuadrantPackage
% Copyright 2013-2016 The University of Queensland.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First, define the particle properties:  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify refractive indices of the medim and particle
n_medium = 1.33; % Water
n_particle = 1.58;


% Particle size
side_length_m =1e-6; % In units of m


% Define the orientations of the cube in (theta, phi) coordinates
Particle_theta=linspace(0,pi/4,3);
Particle_phi = zeros(1,3);



% Next define the measurement setup:  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Objective NA; this defines the trapping beam width. Note, this assumes an
% aberration free objective.
NA=1.25;

% Effective condenser NA
NA_condenser=1.0;

% Polarisation. [ 1 0 ] is plane-polarised along the x-axis, [ 0 1 ] is
% y-polarised, and [ 1 -i ] and [ 1 i ] are circularly polarised.
polarisation = [ 1 -i ];

% Define the spatial mode of the incident light, in the Laguerre-Gauss
% basis. A Gaussian profile is [0, 0].
lg_mode=[0 1]; 


% The vacuum wavelength in m. In this code, we use units of the medium
% wavelength, and used to convert from SI
lambda0=1064e-9;

% To normalize the shot-noise limit, define measured power in W
Power=1e-3; 


% What do we want to calculate:  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Which axis do we measure along: x=1, y=2, z=3
axis=1;


% All needed parameters are now defined.  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% Convert size into medium wavelength units
side_length = side_length_m/lambda0*n_medium;


% Relative index
n_relative = n_particle/n_medium;


% Location of the focal point relative to the particle. These are the
% [ x y z ] coordinates.
beam_offset = [ 0 0 0];


% Specify the beam width. 
beam_angle = asin(NA/n_medium)*180/pi;
w0 = lg_mode_w0( lg_mode, beam_angle );



% To give all measurements in wavelengths in the surrounding medium:
wavelength = 1;
k = 2*pi/wavelength;
k_particle = 2*pi/wavelength*n_relative;


% To what order do we expand the fields?
Nmax = ka2nmax(k*side_length); %nmax for a centered cube with side length radius...
Nmax_medium=Nmax; %these are the same using this method.
Nmax_particle=ka2nmax(k*side_length*n_relative); %this is the internal refractive index in this system



% Create the trapping field:
% Use this for focused Gaussian field
[n,m,a0,b0] = bsc_pointmatch_farfield(Nmax,1,[ lg_mode w0 1 polarisation 90 beam_offset ]);
[a,b,n,m] = make_beam_vector(a0,b0,n,m);


% Normalize total power of wave sum to 1.
pwr = sqrt(sum( abs(a).^2 + abs(b).^2 ));
a=a/pwr;
b=b/pwr;

%********* Insert tmatrix here *********%
T=tmatrix_pm_cube(Nmax,Nmax_medium,Nmax_particle,k,k_particle,side_length);
%*********                     *********%


% Define the grid of points used in measurement
theta0=linspace(0,asin(NA_condenser/n_medium),150);% We only need to include points with theta<theta_max
N_phi=160;% Make this a divisible by 4, so that each quadrant includes the same number of points.
phi0=pi/N_phi:2*pi/N_phi:2*pi;

[theta,phi]= meshgrid(theta0,phi0);
dA=sin(theta(:))*(theta0(2)-theta0(1))*(phi0(2)-phi0(1)); % The area of each grid point


% Calculate the matrices which transfer between expansion coefficients and
% Electric field amplitudes
[A_p,B_q]=farfield_matrix2(Nmax,theta(:),phi(:));


% Define the particle displacements for which we calculate the signal 
dx = linspace(-2.5,2.5,251);


% Generally, the trapping position is not centred along the z axis. We could place the particle at the trap point,
% but instead here place the particle at the beam centre;
zeq=0;


% Now work out spherical coordinates for the displacements:
if(axis==1)
    [rt,theta_dx,phi_dx]=xyz2rtp(dx,0,zeq); % Displacement along (x,y,z) axis
elseif(axis==2)
    [rt,theta_dx,phi_dx]=xyz2rtp(0,dx,zeq); 
elseif(axis==3)
    [rt,theta_dx,phi_dx]=xyz2rtp(0,0,dx);   
    
    if(length(phi_dx)==1) %This seems to be a problem
        phi_dx=phi_dx*ones(size(theta_dx));
    end
else
    warning('axis must be 1, 2, or 3.')
    break
end

% Define these here;
Itot =zeros(length(dx),length(Particle_theta));
QuadX=zeros(length(dx),length(Particle_theta));
QuadY=zeros(length(dx),length(Particle_theta));



for nr = 1:length(dx)
    
    Rot = z_rotation_matrix(theta_dx(nr),phi_dx(nr)); %calculates an appropriate axis rotation off z.
    D = wigner_rotation_matrix(Nmax,Rot);

        
    [A,B] = translate_z(Nmax,rt(nr));
    
        for n_ang = 1:length(Particle_theta)
        % Rotation of the beam to the particle orientation coordinates.
        orientation=z_rotation_matrix(Particle_theta(n_ang),Particle_phi(n_ang));
        D2 = wigner_rotation_matrix(Nmax,orientation);    
        
        a2 = D2*(  D'*A * D*a +  D'*B * D*b ); % Wigner matricies here are hermitian. Therefore in MATLAB the D' operator is the inverse of D.
        b2 = D2*(  D'*A * D*b +  D'*B * D*a ); % In MATLAB operations on vectors are done first, therefore less calculation is done on the matricies.
    
        pq = T * [ a2; b2 ];
        p = pq(1:length(pq)/2);
        q = pq(length(pq)/2+1:end);
    
        % Now translate the scattered field back into the original reference frame
        [A2,B2] = translate_z(Nmax,-rt(nr));
        p2 = D'*(  A2 * D*(D2'*p) +  B2 * D*(D2'*q) );
        q2 = D'*(  A2 * D*(D2'*q) +  B2 * D*(D2'*p) );

        % Calculate the electric field for these coefficients
        E=reshape(A_p*(a+2*p2)+B_q*(b+2*q2),length(theta(:)),2);


        % This is total current and the two subtraction signals.
        Itot(nr,n_ang) =sum(sum(abs(E).^2,2).*dA);
        QuadX(nr,n_ang)=sum(sum(abs(E).^2.*(sign(cos(phi(:)))*ones(1,2)),2).*dA);%/Itot(nr);
        QuadY(nr,n_ang)=sum(sum(abs(E).^2.*(sign(sin(phi(:)))*ones(1,2)),2).*dA);%/Itot(nr);

        end
    
end
% When finding the smallest resolvable displacement, include only the linear region in dx 
Incl=abs(dx)<0.2;
pow=mean(Itot(Incl));

LinfitX=polyfit(dx(Incl)',QuadX(Incl,1)/pow,1);
GX=LinfitX(1);
% LinfitY=polyfit(dx(Incl)',QuadY(Incl,1),1);
% GY=LinfitY(1);
LinfitZ=polyfit(dx(Incl)',Itot(Incl,1)/pow,1);
GZ=LinfitZ(1);

% Minimum resolvable displacement, if a noiseless 100% efficient measurement is performed.
%This is given by 1/(GX*sqrt(n)), with n the measured photon number and GX in units m^-1 

N_photflux=Power*lambda0/(6.63e-34*3.00e8);
dx_min=abs(lambda0/(n_medium*GX *sqrt(N_photflux))); % In units of: m Hz^{-1/2}
dz_min=abs(lambda0/(n_medium*GZ *sqrt(N_photflux)));

% plot(dx*lambda0/n_medium,[QuadX;QuadY;Itot-mean(Itot)]/pow)
plot(dx*lambda0/n_medium,QuadX/pow)
toc