function [A_p,B_q] = farfield_matrix(n,m,theta,phi)
% farfield_matrix.m
% Finds matrices which return the far field from VSWF expansion coefficients
%
% usage
% [A_p,B_q] = farfield_matrix(n,m,theta,phi)
%
% then, the electric field is given by
% E=reshape(A_p*(a+2*p2)+B_q*(b+2*q2),length(theta(:)),3);
%
% each row of E is the field (in spherical coordinates) in the
% (theta,phi) direction (assuming a distance scaling factor of kr)
%
%
%
% This file is part of the Quadrant Detection Package version 1.1. 
% M. A. Taylor and W. P. Bowen, “A computational tool to characterize
% particle tracking measurements in optical tweezers,” J. Opt. 15, 085701
% (2013).  
% https://github.com/michael-a-taylor/QuadrantPackage
% Copyright 2013-2016 The University of Queensland.
%
%
%
% The file is modified from "farfield.m" in the Optical tweezers toolbox 1.2,
% and calls on other functions from that software package.
%
% For details of the Optical tweezers toolbox, see:
% http://www.physics.uq.edu.au/people/nieminen/software.html
%(Copyright 2006-2013 The University of Queensland.)
%

[theta,phi] = matchsize(theta,phi);

A_p = zeros(length(theta)*3,length(n));
B_q = zeros(length(theta)*3,length(n));

for nn = 1:length(n)
   
   [B,C,P] = vsh(n(nn),m(nn),theta,phi);
   Nn = 1/sqrt(n(nn)*(n(nn)+1));
      
   
   A_p(:,nn) =  Nn * (-i)^(n(nn)+1)*reshape(C,1,length(theta)*3);
   B_q(:,nn) =  Nn * (-i)^n(nn)*reshape(B,1,length(theta)*3);      
   
   
end

return