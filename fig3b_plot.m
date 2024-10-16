clear all

DeltaA = 6.5e-18; % m^2
N = 1000;
k = 1; % 1/s
a = 1e-5; % m
eta = 1e-3; % Pa*s
beta = 2.5e20; % 1/J
D = 1.5e-10; % m^2/s
K = 5e-14; % m^2

vvec = linspace(1e-9,2.5e-6,1e3); % m/s
Gvec = logspace(-3,2,1e3);
v = vvec'*ones(1,length(Gvec));
G = ones(length(vvec),1)*Gvec;
v_ = vvec*1e6; % um/s

% precision ratio, P_chem^exo/P_mech
q = beta*eta*v*a^2*DeltaA/K;
R = G*sqrt(D/a^2)./q/sqrt(N*k);

% bound color range
Rmax = 1e3;
Rmin = 1e-3;
R(R>Rmax) = Rmax;
R(R<Rmin) = Rmin;

% create theory region for G <= 0 (R <= 0)
G2vec = linspace(-.2,0,10);
R2 = Rmin*ones(length(v_),length(G2vec));

% create colormap
Z = 1e2;
cmin = log10(Rmin);
cmax = log10(Rmax);
cmap1 = [linspace(0,1,Z)' linspace(0,1,Z)' ones(Z,1)]; % blue to white
cmap2 = [ones(Z,1) linspace(1,0,Z)' linspace(1,0,Z)']; % white to red
cmap = [cmap1;cmap2];

fs = 30; fs2 = 26;
ms = 20;
lw = 2;
w = .03; % um/s (error bar cap width)

% theory
figure(1); clf
subplot(10,8,[2:8 10:16 18:24 26:32 34:40 42:48 50:56 58:64])
imagesc(v_,Gvec,log10(R)')
set(gca,'ydir','normal','yscale','log','fontsize',fs2,...
    'xticklabels',[])
xlim([0 max(v_)])
ylim([min(Gvec) max(Gvec)])
ylabel('Gradient, $G = g\sqrt{a^5/c}\quad\,\,$',...
    'interpreter','latex','fontsize',fs)
title('${\cal P}_{\rm chem}^{\rm exo}/{\cal P}_{\rm mech}$',...
    'interpreter','latex','fontsize',28)
caxis([cmin cmax])
colormap(cmap)

subplot(10,8,66:72)
imagesc(v_,G2vec,log10(R2)')
set(gca,'ydir','normal','fontsize',fs2,'ytick',[-.2 0])
xlim([0 max(v_)])
ylim([min(G2vec) max(G2vec)])
xlabel('Fluid velocity, $v$ ($\mu$m/s)','interpreter','latex',...
    'fontsize',fs)
caxis([cmin cmax])
colormap(cmap)

% Moon 2023 data
load('fig3b.mat')
% flow only -- upstream
v0 = 1.5; % um/s
dv0 = .05; % um/s
G0bar = 0;
% chemical only -- downstream
v1 = 0; % um/s
% chemical vs flow, small-G region -- upstream
v2 = 1.5; % um/s
dv2 = .05; % um/s
% chemical vs flow, large-G region -- downstream
v3 = 1.5; % um/s
dv3 = .05; % um/s
% chemical with flow, small-G region (saturation) -- nowhere
v4 = 1.5; % um/s
dv4 = .05; % um/s

subplot(10,8,[2:8 10:16 18:24 26:32 34:40 42:48 50:56 58:64])
hold on
plot([v1 v1],G1bar+dG1*[-1 1],'k-','linewidth',lw)
plot([v1-w v1+w],(G1bar-dG1)*[1 1],'k-','linewidth',lw)
plot([v1-w v1+w],(G1bar+dG1)*[1 1],'k-','linewidth',lw)
plot([v2 v2],G2bar+dG2*[-1 1],'k-','linewidth',lw)
plot([v2-w v2+w],(G2bar-dG2)*[1 1],'k-','linewidth',lw)
plot([v2-w v2+w],(G2bar+dG2)*[1 1],'k-','linewidth',lw)
plot([v3 v3],G3bar+dG3*[-1 1],'k-','linewidth',lw)
plot([v3-w v3+w],(G3bar-dG3)*[1 1],'k-','linewidth',lw)
plot([v3-w v3+w],(G3bar+dG3)*[1 1],'k-','linewidth',lw)
plot(v1,G1bar,'ko','markersize',ms,'linewidth',lw,'markerfacecolor','r')
plot(v2,G2bar,'ko','markersize',ms,'linewidth',lw,'markerfacecolor','b')
plot(v3,G3bar,'ko','markersize',ms,'linewidth',lw,'markerfacecolor','r')

subplot(10,8,66:72)
hold on
plot(v0,G0bar,'ko','markersize',ms,'linewidth',lw,'markerfacecolor','b')
box on

plot(v4,G4bar,'ko','markersize',ms,'linewidth',lw,'markerfacecolor','w')
box on
print(gcf,'-depsc','fig3b.eps')

% for colorbar
figure(2); clf
imagesc(v_,Gvec,log10(R)')
set(gca,'ydir','normal','yscale','log','fontsize',fs2,...
    'xticklabels',[])
xlim([0 max(v_)])
ylim([min(Gvec) max(Gvec)])
caxis([cmin cmax])
colormap(cmap)
ch = colorbar;
ch.TickLabels = {'≤10^{-3}','10^{-2}','10^{-1}','10^0',...
    '10^1','10^2','≥10^3'};
print(gcf,'-depsc','fig3b_colorbar.eps')
