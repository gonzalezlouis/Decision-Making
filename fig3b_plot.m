clear all

DeltaA = 6.5e-18; % m^2
N = 1000;
k = 1; % 1/s
a = 1e-5; % m
eta = 1e-3; % Pa*s
beta = 2.5e20; % 1/J
D = 1.5e-10; % m^2/s
K = 5e-14; % m^2

vvec = linspace(1e-9,3e-6,1e3); % m/s
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
subplot(30,8,[2:8 10:16 18:24 26:32 34:40 42:48 50:56 58:64 ...
    66:72 74:80 82:88 90:96 98:104 106:112 114:120 122:128 ...
    130:136 138:144 146:152 154:160 162:168 170:176 178:184])
imagesc(v_,Gvec,log10(R)')
set(gca,'ydir','normal','yscale','log','fontsize',fs2,...
    'xticklabels',[],'ytick',10.^(-2:2:2))
xlim([0 max(v_)])
ylim([min(Gvec) max(Gvec)])
ylabel('Gradient, $G = g\sqrt{a^5/c}\quad\quad\,$',...
    'interpreter','latex','fontsize',fs)
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
% chemical only -- up-gradient
v1 = 0; % um/s
% chemical vs flow, small-G region -- upstream
v2 = 1.5; % um/s
dv2 = .05; % um/s
% chemical vs flow, large-G region -- up-gradient
v3 = 1.5; % um/s
dv3 = .05; % um/s
% chemical with flow, small-G region (saturation) -- neither
v4 = 1.5; % um/s
dv4 = .05; % um/s

subplot(30,8,[2:8 10:16 18:24 26:32 34:40 42:48 50:56 58:64 ...
    66:72 74:80 82:88 90:96 98:104 106:112 114:120 122:128 ...
    130:136 138:144 146:152 154:160 162:168 170:176 178:184])
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
plot(v4,G4bar,'ko','markersize',ms,'linewidth',lw,'markerfacecolor','w')

% Our data
L = 1e3; % um
x = linspace(0,L,1e3); % um
dx = x(2)-x(1);
d = 100; % um -- excluded edge width
[ig,i1] = min(abs(x-d));
[ig,i2] = min(abs(x-(L-d)));

% chemical only -- up-gradient
v5 = 0;
c1 = 0; % 1/um^3
c2 = 5*.6; % 1/um^3
c = c1 + x*(c2-c1)/L; % 1/um^3
g = (c2-c1)/L; % 1/um^4
G = g*sqrt(a^5./c);
G5bar = sum(G(i1:i2))/(L-2*d);
dG5 = sqrt(sum((G(i1:i2)).^2)/(L-2*d) - G5bar^2);

% flow only -- upstream
v6 = 2.5; % um/s
G6bar = 0;

% chemical vs flow, small G -- upstream
v7 = 2.5; % um/s
c1 = 2.5*.6; % 1/um^3
c2 = 5*.6; % 1/um^3
Du = D*(1e6)^2; % um^2/s
c = c1 + (c2-c1)*((exp(v7*x/Du)-1)/(exp(v7*L/Du)-1)); % 1/um^3
g = v7*(c2-c1)*exp(v7*x/Du)/Du/(exp(v7*L/Du)-1); % 1/um^4
G = g.*sqrt(a^5./c);
G7bar = sum(G(i1:i2))/(L-2*d);
dG7 = sqrt(sum((G(i1:i2)).^2)/(L-2*d) - G7bar^2);

% chemical vs flow, large G -- up-gradient
v8 = 2.5; % um/s
c1 = 0; % 1/um^3
c2 = 10*.6; % 1/um^3
Du = D*(1e6)^2; % um^2/s
c = c1 + (c2-c1)*((exp(v8*x/Du)-1)/(exp(v8*L/Du)-1)); % 1/um^3
g = v8*(c2-c1)*exp(v8*x/Du)/Du/(exp(v8*L/Du)-1); % 1/um^4
G = g.*sqrt(a^5./c);
G8bar = sum(G(i1:i2))/(L-2*d);
dG8 = sqrt(sum((G(i1:i2)).^2)/(L-2*d) - G8bar^2);

subplot(30,8,[2:8 10:16 18:24 26:32 34:40 42:48 50:56 58:64 ...
    66:72 74:80 82:88 90:96 98:104 106:112 114:120 122:128 ...
    130:136 138:144 146:152 154:160 162:168 170:176 178:184])
plot([v5 v5],G5bar+dG5*[-1 1],'k-','linewidth',lw)
plot([v5-w v5+w],(G5bar-dG5)*[1 1],'k-','linewidth',lw)
plot([v5-w v5+w],(G5bar+dG5)*[1 1],'k-','linewidth',lw)
plot(v5,G5bar,'k^','markersize',ms,'linewidth',lw,'markerfacecolor','r')
plot([v8 v8],[G8bar*.5 G8bar+dG8],'k-','linewidth',lw) % down-arrow pt 1
plot([v8 v8+w],[G8bar*.5 G8bar*.6],'k-','linewidth',lw) % down-arrow pt 2
plot([v8 v8-w],[G8bar*.5 G8bar*.6],'k-','linewidth',lw) % down-arrow pt 3
plot([v8-w v8+w],(G8bar-dG8)*[1 1],'k-','linewidth',lw)
plot([v8-w v8+w],(G8bar+dG8)*[1 1],'k-','linewidth',lw)
plot(v8,G8bar,'k^','markersize',ms,'linewidth',lw,'markerfacecolor','r')
plot([v7 v7],[G7bar*.5 G7bar+dG7],'k-','linewidth',lw) % down-arrow pt 1
plot([v7 v7+w],[G7bar*.5 G7bar*.6],'k-','linewidth',lw) % down-arrow pt 2
plot([v7 v7-w],[G7bar*.5 G7bar*.6],'k-','linewidth',lw) % down-arrow pt 3
plot([v7-w v7+w],(G7bar-dG7)*[1 1],'k-','linewidth',lw)
plot([v7-w v7+w],(G7bar+dG7)*[1 1],'k-','linewidth',lw)
plot(v7,G7bar,'k^','markersize',ms,'linewidth',lw,'markerfacecolor','b')
box on

subplot(10,8,66:72)
plot(v6,G6bar,'k^','markersize',ms,'linewidth',lw,'markerfacecolor','b')
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
