clear all

DeltaA = 6.5e-18; % m^2
N = 1000;
k = 1; % 1/s
nu = 1; % 1/s
a = 1e-5; % m
eta = 1e-3; % Pa*s
beta = 2.5e20; % 1/J
D = 1.5e-10; % m^2/s
L = 3e-3; % m
K = 1e-13; % m^2

vvec = logspace(-8,-4,1e3); % m/s
rhovec = logspace(9,13,1e3); % 1/m^3
v = vvec'*ones(1,length(rhovec));
rho = ones(length(vvec),1)*rhovec;
v_ = vvec*1e6; % um/s
rho_ = rhovec/1e9; % 1/mm^3

% precision ratio, P_chem^auto/P_mech
epsilon = v*a/D;
q = beta*eta*v*a^2*DeltaA/K;
Pc = epsilon.*sqrt(nu*epsilon./(epsilon+a^2*L*rho));
Pm = q*sqrt(N*k);
R = Pc./Pm;

% create colormap
Z = 1e2;
cmin = min(R(:));
cmax = 1 + 1-cmin;
cmap1 = [linspace(0,1,Z)' linspace(0,1,Z)' ones(Z,1)]; % blue to white
cmap2 = [ones(Z,1) linspace(1,0,Z)' linspace(1,0,Z)']; % white to red
cmap = [cmap1;cmap2];

fs = 30; fs2 = 26;
ms = 20;
lw = 2;

% theory
figure(1); clf
imagesc(v_,rho_,R')
set(gca,'ydir','normal','yscale','log','xscale','log',...
    'fontsize',fs2)
xlim([min(v_) max(v_)])
ylim([min(rho_) max(rho_)])
xlabel('Fluid velocity, $v$ ($\mu$m/s)','interpreter','latex',...
    'fontsize',fs)
ylabel('Cell density, $\rho$ (mm$^{-3}$)','interpreter','latex',...
    'fontsize',fs)
title('${\cal P}_{\rm chem}^{\rm auto}/{\cal P}_{\rm mech}$',...
    'interpreter','latex','fontsize',fs)
caxis([cmin cmax])
colormap(cmap)
colorbar

% Polacheck 2011 data
% low density -- downstream
v1 = [.3 3]; % um/s
rho1 = [50 50]; % 1/mm^3
% high density -- upstream
v2 = [.3 3]; % um/s
rho2 = [250 250]; % 1/mm^3
% Polacheck 2014 data
% higher density, higher velocity -- upstream
v3 = 4.6; % um/2
rho3 = 600; % 1/mm^3

% add 2011 data
hold on
plot(v1,rho1,'ko','markersize',ms,'linewidth',lw,'markerfacecolor','r')
plot(v2,rho2,'ko','markersize',ms,'linewidth',lw,'markerfacecolor','b')

% add 2014 data
plot(v3,rho3,'ks','markersize',ms,'linewidth',lw,'markerfacecolor','b')
box on
print(gcf,'-depsc','fig3a.eps')
