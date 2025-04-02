clear all

a = 1e-5; % m
D = 1.5e-10; % m^2/s

D2c2vec = logspace(-2,4,1e3); % 1/s
Gvec = linspace(0,2,1e3);
D2c2 = D2c2vec'*ones(1,length(Gvec));
G = ones(length(D2c2vec),1)*Gvec;

% precision ratio, P_chem^exo/P_CIL
R = G*sqrt(D)/a./sqrt(D2c2);

% bound color range
Rmax = 1e2;
Rmin = 1e-2;
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
ms = 20;
lw = 2;
w = .03; % (error bar cap width)

% plot theory
figure(1); clf
imagesc(D2c2vec,Gvec,log10(R)')
set(gca,'ydir','normal','xscale','log','fontsize',fs2,...
    'xtick',10.^(-2:2:4))
xlim([min(D2c2vec) max(D2c2vec)])
ylim([min(Gvec) max(Gvec)])
xlabel('Ephrin turnover, $D_2c_2$ (s$^{-1}$)',...
    'interpreter','latex','fontsize',fs)
ylabel('Gradient, $G = g\sqrt{a^5/c}$',...
    'interpreter','latex','fontsize',fs)
title('${\cal P}_{\rm chem}^{\rm exo}/{\cal P}_{\rm CIL}$',...
    'interpreter','latex','fontsize',fs)
caxis([cmin cmax])
colormap(cmap)
ch = colorbar;
ch.Ticks = -2:2;
ch.TickLabels = {'≤10^{-2}','10^{-1}','10^0',...
    '10^1','≥10^2'};

% Lin 2015 data
c = [4.1 8.3 12.4]*.6; % 1/um^3
Deltac = 3.3*.6; % 1/um^3
f = .4;
d = 250; % um
g = f*Deltac/d; % 1/um^4
a = 10; % um
G = g*sqrt(a^5./c);
G1 = G(1);
G2 = G(2);
G3 = G(3);
Dc = 600; % 1/s
Dc_range = [200 2000]; % 1/s

% add data
hold on
plot(Dc_range,[G1 G1],'k-','linewidth',lw)
plot(Dc_range(1)*[1 1],[G1-w G1+w],'k-','linewidth',lw)
plot(Dc_range(2)*[1 1],[G1-w G1+w],'k-','linewidth',lw)
plot(Dc_range,[G2 G2],'k-','linewidth',lw)
plot(Dc_range(1)*[1 1],[G2-w G2+w],'k-','linewidth',lw)
plot(Dc_range(2)*[1 1],[G2-w G2+w],'k-','linewidth',lw)
plot(Dc_range,[G3 G3],'k-','linewidth',lw)
plot(Dc_range(1)*[1 1],[G3-w G3+w],'k-','linewidth',lw)
plot(Dc_range(2)*[1 1],[G3-w G3+w],'k-','linewidth',lw)
plot(Dc,G1,'ko','markersize',ms,'linewidth',lw,'markerfacecolor','r')
plot(Dc,G2,'ko','markersize',ms,'linewidth',lw,'markerfacecolor','r')
plot(Dc,G3,'ko','markersize',ms,'linewidth',lw,'markerfacecolor','b')
box on
print(gcf,'-depsc','fig3c.eps')

