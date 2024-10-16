clear all
% dimensionless sensitivity parameter: G = g*sqrt(a^5/c)

% data from Moon et al, 2023
moon2D = load('moon2D.txt');
moon2E = load('moon2E.txt');
moon2F = load('moon2F.txt');
a = 10; % um

% chemical only (just points in Fig 2D)
x1 = moon2D(9:-1:3,1); % um
c1 = moon2D(9:-1:3,2)*.6; % 1/um^3
g1 = diff(c1)./diff(x1); % 1/um^4
c1_ = (c1(2:end) + c1(1:end-1))/2;
G1 = g1.*sqrt(a^5./c1_);
G1bar = mean(G1);
dG1 = std(G1);

% chemical and flow, large G (just points in right region of Fig 2F)
x3 = moon2F(9:-1:4,1); % um
c3 = moon2F(9:-1:4,2)*.6; % 1/um^3
g3 = diff(c3)./diff(x3); % 1/um^4
c3_ = (c3(2:end) + c3(1:end-1))/2;
G3 = g3.*sqrt(a^5./c3_);
G3bar = mean(G3);
dG3 = std(G3);

% chemical and counter-flow (just points in right region of Fig 2E)
x4 = moon2E(7:13,1); % um
c4 = moon2E(7:13,2)*.6; % 1/um^3
g4 = diff(c4)./diff(x4); % 1/um^4
c4_ = (c4(2:end) + c4(1:end-1))/2;
G4 = g4.*sqrt(a^5./c4_);
G4bar = -mean(G4); % negative because gradient points upstream here
dG4 = std(G4);

% chemical and flow, small G: need to use raw data from Fig 2F
x = load('moon2F_raw_x.txt');
n = load('moon2F_raw_y.txt');
n = n - min(n(:)); % subtract background (lowest n value must be zero)
[x,j] = sort(x);
n = n(j,:);
nbar = mean(n,2);

% left region used in Fig 2F
x1 = 125;
x2 = 420;
[ig,j1] = min(abs(x-x1));
[ig,j2] = min(abs(x-x2));
x_ = x(j1:j2);
cs = n(j1:j2,:)*10*.6; % nM

a = 10; % um
for i = 1:6
    c = cs(:,i);
    g = diff(c)./diff(x_); % 1/um^4
    cmid = (c(2:end) + c(1:end-1))/2;
    G(:,i) = g.*sqrt(a^5./cmid);
end

% avg over whole region (sd from multiple trials)
Gs = mean(G,1);
G2bar = mean(Gs);
dG2 = std(Gs);

figure(1); clf
h = plot(x,n,'.',x,nbar,'k-',[-400 1400],[0 0],'k-',...
    [x1 x1],[-.2 1.2],'r-',[x2 x2],[-.2 1.2],'r-');
set(h(7),'linewidth',2)
xlabel('x (\mum)')
ylabel('normalized concentration')

save('fig3b.mat')