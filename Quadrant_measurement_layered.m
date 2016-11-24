tic

% This file will calculate the measurement signal and minimum resolvable
% displacement for a *layered* spherical particle in optical tweezers, with
% particle tracking via a quadrant detector at the back-focal plane of a
% condenser. 

% This file is part of the Quadrant Detection Package version 1.1. 
% M. A. Taylor and W. P. Bowen, “A computational tool to characterize
% particle tracking measurements in optical tweezers,” J. Opt. 15, 085701
% (2013).  
% https://github.com/michael-a-taylor/QuadrantPackage
% Copyright 2013-2016 The University of Queensland.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First, define the particle properties:  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify refractive indices of the medim, and the particle core and shell
n_medium = 1.33; % Water

% Core material
n_core = 2.3;

% Shell refractive index. To use multiple layers, the n_shell variable
% should be a row vector containing all the refractive indices of the shells.
% n_shell = [1.4, 1.6, 1.4];
n_shell = 1.78;

% Particle core size. Specify this in m
radius_core = 250e-9; 


% Define the shell thicknesses in m.
% shell_thickness_m = [25e-9, 25e-9, 25e-9]; 
shell_thickness = 230e-9;


% Next define the measurement setup:  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Objective NA; this defines the trapping beam width. Note, this assumes an
% aberration free objective.
NA=1.25;

% Effective condenser NA
NA_condenser=1.25;

% Polarisation. [ 1 0 ] is plane-polarised along the x-axis, [ 0 1 ] is
% y-polarised, and [ 1 -i ] and [ 1 i ] are circularly polarised.
polarisation = [ 1 0 ];


% Define the spatial mode of the incident light, in the Laguerre-Gauss
% basis. A Gaussian profile is [0, 0].
lg_mode=[0 0]; 


% The vacuum wavelength in m. In this code, we use units of the medium
% wavelength, and used to convert from SI
lambda0=1064e-9;


% To normalize the shot-noise limit, define measured power in W
Power=1e-3; 



% What do we want to calculate:  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Which axis do we measure along: x=1, y=2, z=3
axis=1;

% Is the particle axially centered at the stable trap point (=1) or the beam focus (=0)?
Trap=1;

% All needed parameters are now defined.  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% Creating the trapping field  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Location of the focal point relative to the particle. These are the
% [ x y z ] coordinates.
beam_offset = [ 0 0 0];


% Specify the beam width. 
beam_angle = asin(NA/n_medium)*180/pi;
w0 = lg_mode_w0( lg_mode, beam_angle );



% Convert sizes into medium wavelength units
radius = radius_core/lambda0*n_medium;
shell_thick = shell_thickness/lambda0*n_medium;


% To give all measurements in wavelengths in the surrounding medium:
wavelength = 1;
k = 2*pi/wavelength;

% To what order do we expand the fields?
Nmax = ka2nmax(k*(radius+sum(shell_thickness)));
    if Nmax < 12
        Nmax = 12;
    end


% Create the trapping field:
% Use this for focused Gaussian field
[n,m,a0,b0] = bsc_pointmatch_farfield(Nmax,1,[ lg_mode w0 1 polarisation 90 beam_offset ]);
[a,b,n,m] = make_beam_vector(a0,b0,n,m);


% Normalize total power of wave sum to 1.
pwr = sqrt(sum( abs(a).^2 + abs(b).^2 ));
a=a/pwr;
b=b/pwr;

%********* Insert tmatrix here *********%
T=tmatrix_mie_layered(Nmax,k,k*[n_core,n_shell]/n_medium,[radius,radius+shell_thick]);
%*********                     *********%


if(Trap)
    % Calculate axial trapping point %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Generally, the trapping position is not centred along the z axis. 

    % Specify points at which to evaluate the force.
    % Note that these lengths are in units of the medium wavelength.
    z = linspace(-5,4,100);
    fz = zeros(size(z));

    %calculate the force along z
    for nz = 1:length(z)
    
        [A,B] = translate_z(Nmax,z(nz));
        a2 = ( A*a + B*b );
        b2 = ( A*b + B*a );
    
        pq = T * [ a2; b2 ];
        p = pq(1:length(pq)/2);
        q = pq(length(pq)/2+1:end);
    
        fz(nz) = force_z(n,m,a2,b2,p,q);
    
    end

    % Locate the trapping point
    maxforce=z(fz==max(fz));
    zeroindex=find(fz<0&z>maxforce,1);

    if length(zeroindex)~=0
        % fit to third order polynomial the local points. (only works when
        % dz sufficiently small)
        pz=polyfit(z(max([zeroindex-2,1]):min([zeroindex+2,length(z)])),fz(max([zeroindex-2,1]):min([zeroindex+2,length(z)])),3);
        root_z=roots(pz); %find roots of 3rd order poly.
    
        dpz=[3*pz(1),2*pz(2),1*pz(3)]; %derivative of 3rd order poly.
    
        real_z=root_z(imag(root_z)==0); % finds real roots only.
    
        rootsofsign=polyval(dpz,real_z); %roots that are stable
        zeq=real_z(rootsofsign<0); %there is at most 1 stable root. critical roots give error.
        try
            zeq=zeq(abs(zeq-z(zeroindex))==min(abs(zeq-z(zeroindex))));
        end
    else
        zeq=[];
    end

    if length(zeq)==0
        warning('No axial equilibrium in range!')
        zeq=0;
    end


else
    % Alternatively, we could simply place the particle at the beam centre;
    zeq=0;
end






% Calculate the detection signal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the grid of points used in measurement
theta0=linspace(0,asin(NA_condenser/n_medium),150);% We only need to include points with theta<theta_max
N_phi=160;% Make this a divisible by 4, so that each quadrant includes the same number of points.
phi0=pi/N_phi:2*pi/N_phi:2*pi;

[theta,phi]= meshgrid(theta0,phi0);
dA=sin(theta(:))*(theta0(2)-theta0(1))*(phi0(2)-phi0(1)); % The area of each grid point

% Calculate the matrices which transfer from expansion coefficients to
% Electric field amplitudes
[A_p,B_q]=farfield_matrix2(Nmax,theta(:),phi(:));


% Define the particle displacements for which we calculate the signal, in
% units of the medium wavelength
dx = linspace(-2,2,201); 



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
Itot =zeros(size(dx));
QuadX=zeros(size(dx));
QuadY=zeros(size(dx));



for nr = 1:length(dx)
    
    Rot = z_rotation_matrix(theta_dx(nr),phi_dx(nr)); %calculates an appropriate axis rotation off z.
    D = wigner_rotation_matrix(Nmax,Rot);
    
    [A,B] = translate_z(Nmax,rt(nr));
    a2 = D'*(  A * D*a +  B * D*b ); % Wigner matricies here are hermitian. Therefore in MATLAB the D' operator is the inverse of D.
    b2 = D'*(  A * D*b +  B * D*a ); % In MATLAB operations on vectors are done first, therefore less calculation is done on the matricies.
    
    pq = T * [ a2; b2 ];
    p = pq(1:length(pq)/2);
    q = pq(length(pq)/2+1:end);
    
    % Now translate the scattered field back into the original reference frame
    [A2,B2] = translate_z(Nmax,-rt(nr));
    p2 = D'*(  A2 * D*p +  B2 * D*q );
    q2 = D'*(  A2 * D*q +  B2 * D*p );

    % Calculate the electric field for these coefficients
    E=reshape(A_p*(a+2*p2)+B_q*(b+2*q2),length(theta(:)),2);


    % This is total current and the two subtraction signals.
    Itot(nr) =sum(sum(abs(E).^2,2).*dA);
    QuadX(nr)=sum(sum(abs(E).^2.*(sign(cos(phi(:)))*ones(1,2)),2).*dA);%/Itot(nr);
    QuadY(nr)=sum(sum(abs(E).^2.*(sign(sin(phi(:)))*ones(1,2)),2).*dA);%/Itot(nr);

end

% When finding the smallest resolvable displacement, include only the linear region in dx 
Incl=abs(dx)<0.2;
pow=mean(Itot(Incl));

LinfitX=polyfit(dx(Incl),QuadX(Incl)/pow,1);
GX=LinfitX(1);
% LinfitY=polyfit(dx(Incl),QuadY(Incl),1);
% GY=LinfitY(1);
LinfitZ=polyfit(dx(Incl),Itot(Incl)/pow,1);
GZ=LinfitZ(1);

% Minimum resolvable displacement, if a noiseless 100% efficient measurement is performed.
%This is given by 1/(GX*sqrt(n)), with n the measured photon number and GX in units m^-1 

N_photflux=Power*lambda0/(6.63e-34*3.00e8);
dx_min=abs(lambda0/(n_medium*GX *sqrt(N_photflux))); % In units of: m Hz^{-1/2}
dz_min=abs(lambda0/(n_medium*GZ *sqrt(N_photflux)));

plot(dx*lambda0/n_medium,[QuadX;QuadY;Itot-mean(Itot)]/pow)
xlim([-1.6,1.6]*1e-6)
toc