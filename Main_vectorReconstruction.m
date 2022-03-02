%Test the reconstruction result when a circular domain is used both for the
%forward and inverse  problem
%The available data is both the longitudinal and transverse integral
%and cimmino algorithm is used for the estimation of the electric field
%please add regularization toolbox: http://www.imm.dtu.dk/~pcha/Regutools/ to the matlab path


clear all; close all; clc

%load forward and inverse triangular mesh
load meshes %(H,g,node)



%% Forward Problem: field and potentials
%electric field from the FEM solution of potential distribution
%estimate the gradient of the potentials (\nabla u)

%load dipole source
load ElectricDipole

%Estimate potential on the domain and Electric field
%load lead field matrix for the estimation of the potetials on the domain
load LeadFieldMatrix

Pot = LF_frw*dip;
Grad = GradientEst(Node_frw,Pot); 
E_frw   = -Grad;                        %E=-\nabla u

M = InterpolateMatrix2d(H_frw,g_frw,g_inv);

PE_frw = M*E_frw;                          %Project field to the inverse mesh

%% load projection matrices

load TransverseRay
load LongitudinalRay


%% Measurements - Forward projects
b = Ray_frw*E_frw(:);  % Longitudinal measurements
c = TRay_frw*E_frw(:); % Transverse measurements

%Observations
ob =[b; c];

%% Inverse Solution
%projection matrix
A = [Ray_inv;TRay_inv];
%%Estimate solution using ART
E_hat = cimmino(A,ob,15);

%Results
figure;
set(gcf, 'Units','centimeters', 'Position',[5 5 14 13])

% -- compare and plot results
axes('position',[0.46 0.975 0.10 0.45])
texthvc(0,0,'Magnitude of the field', 'cl',[0 0 0])
axis off
axes('position',[0.1 0.57 0.32 0.34]);
Plotinvsolnode(abs(PE_frw(:,1)+1i*PE_frw(:,2)),g_inv,H_inv);
title('Test Case')
axis off
axis equal
axes('position',[0.15 0.50 0.25 0.08]);
colorbar('location','South')
caxis([0 max(abs(PE_frw(:,1)+1i*PE_frw(:,2)))])  
axis off


E_hat = reshape(E_hat,length(E_hat)/2,2);

axes('position',[0.52 0.57 0.32 0.34])
Plotinvsolnode(abs(E_hat(:,1)+1i*E_hat(:,2)),g_inv,H_inv);
caxis([0 max(abs(E_hat(:,1)+1i*E_hat(:,2)))])  
title(['Reconstruction']); 
axis equal
axis off

axes('position',[0.57 0.50 0.25 0.08]);
colorbar('location','South')
caxis([0 max(abs(E_hat(:,1)+1i*E_hat(:,2)))])  
axis off



axes('position',[0.46 0.46 0.10 0.45])
texthvc(0,0,'(Unit-length) field lines', 'cl',[0 0 0])
axis off

axes('position',[0.1 0.09 0.32 0.34]);
PnEfrw = NormalizeField(PE_frw);
quiver(g_inv(:,1),g_inv(:,2),PnEfrw(:,1),PnEfrw(:,2));
axis equal
axis off
axes('position',[0.52 0.09 0.32 0.34])
nE_hat = NormalizeField([E_hat(:,1) E_hat(:,2)]);
quiver(g_inv(:,1),g_inv(:,2),nE_hat(:,1),nE_hat(:,2));
axis equal
axis off


