clear all

DeltaA = 6.5e-18; % m^2
N = 1000;
k = 1; % 1/s
nu = 6; % 1/s
a = 1e-5; % m
eta = 1e-3; % Pa*s
beta = 2.5e20; % 1/J
D = 1.5e-10; % m^2/s
L = 3e-3; % m
K = 1e-13; % m^2

vvec = logspace(-7,-6,1e3); % m/s
rhovec = logspace(11,12,1e3); % 1/m^3
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

% bound color range
Rmax = 10^.5;
Rmin = 10^(-.5);
R(R>Rmax) = Rmax;
R(R<Rmin) = Rmin;

% create colormap
Z = 1e2;
cmin = log10(Rmin);
cmax = log10(Rmax);
cmap1 = [linspace(0,1,Z)' linspace(0,1,Z)' ones(Z,1)]; % blue to white
cmap2 = [ones(Z,1) linspace(1,0,Z)' linspace(1,0,Z)']; % white to red
cmap = [cmap1;cmap2];

fs = 30; fs2 = 26;
ms = 20*2.5;
lw = 2*2.5;
w = 30; % (error bar cap width)

% theory
figure(1); clf
imagesc(v_,rho_,log10(R)')
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
ch = colorbar;
ch.Ticks = -.5:.5:.5;
ch.TickLabels = {'≤10^{-0.5}','10^0','≥10^{0.5}'};
set(gca,'ydir','normal','yscale','log','xscale','log',...
    'fontsize',fs2,'xtick',10.^[-1 0])

% Munson 2013 data
v1 = .7; % um/s
rho1 = 300; % 1/mm^3
v2 = .5; % um/s
v2_range = [.2 .8]; % um/s
rho2 = 300; % 1/mm^3

% add data
hold on
plot(v1,rho1,'kd','markersize',ms,'linewidth',lw,'markerfacecolor','r')
plot(v2_range,[rho2 rho2],'k-','linewidth',lw)
plot(v2_range(1)*[1 1],[rho2-w rho2+w],'k-','linewidth',lw)
plot(v2_range(2)*[1 1],[rho2-w rho2+w],'k-','linewidth',lw)
plot(v2,rho2,'kd','markersize',ms,'linewidth',lw,'markerfacecolor','r')
box on

print(gcf,'-depsc','fig3a_inset.eps')
