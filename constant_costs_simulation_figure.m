% This code simulates the evolution of the rate of onset of resistance,
% zeta.

%% Plot 1 (host mortality costs and parasite mortality virulence)

% Define parameters:
version=2;
a0=5;
c1a=0.1;
c2a=2;
b0=0.1;
c1b=0.1;
c2b=2;
f=1;
h=1;
beta=1;
alpha=1;
delta=0.25;
q=1;
gamma=0;

t_max=100;
zetamin=0;
zetamax=1;
zetastart=0.1;
res0=51;
nevol=4000;

% Set up vector to use later:
ZETA=[];

% Set up initial conditions for a monomorphic population:
strain_total = 1;
init_pop=[0.1,0.1,0];
Zeta = linspace(zetamin,zetamax,res0);
initial = find(Zeta>=zetastart,1);
zeta_start = Zeta(initial);
index_start = initial;

% Allow zeta to evolve:
[~,~,~,~,ZETAnew,~,~,~] = simulation_function(t_max,a0,c1a,c2a,b0,c1b,c2b,beta,alpha,delta,zetamin,zetamax,zeta_start,h,f,q,gamma,init_pop,strain_total,index_start,res0,nevol,version);
ZETA=[ZETA;ZETAnew];

% Plot the evolutionary trajectory:
subplot(2,1,1)
aa0=1;
ZETA0=log10(ZETA);
ZETA0(ZETA0<-aa0)=-aa0;
ZETA0=(ZETA0+aa0)/aa0;
imagesc(ZETA0);set(gca,'ydir','normal')
map=colormap('gray');
map=flipud(map);
colormap(map);
set(gca,'xtick',1:(res0-1)/2:res0,'xticklabel',round(Zeta(1:(res0-1)/2:res0)*10000)/10000,'fontsize',12);
set(gca,'ytick',[0,1000,2000,3000,4000]);
ylabel('Evolutionary time','interpreter','latex','fontsize',16)
xlabel('Rate of onset of resistance, $\zeta$','interpreter','latex','fontsize',16)
pbaspect([1,1.5,1])
hold on

% Plot the analytically determined value of the singular strategy:
singstrat=85/192;
singstratval=1+(singstrat/(1/(res0-1)));
singstrat_x=singstratval+zeros(4000,1);
singstrat_y=linspace(1,4000,4000);
plot(singstrat_x,singstrat_y,'--r','linewidth',2)
text(1,3900,'A','fontsize',16)

%% Plot 2 (host fecundity costs and parasite sterility virulence)

% Define parameters:
version=1;
a0=5;
c1a=0.1;
c2a=2;
b0=0.1;
c1b=0.1;
c2b=2;
f=0.5;
h=0.75;
beta=1;
alpha=0;
delta=0;
q=1;
gamma=0;

t_max=100;
zetamin=0;
zetamax=1;
zetastart=0.1;
res0=51;
nevol=4000;

% Set up vector to use later:
ZETA=[];

% Set up initial conditions for a monomorphic population:
strain_total = 1;
init_pop=[0.1,0.1,0];
Zeta = linspace(zetamin,zetamax,res0);
initial = find(Zeta>=zetastart,1);
zeta_start = Zeta(initial);
index_start = initial;

% Allow zeta to evolve:
[~,~,~,~,ZETAnew,~,~,~] = simulation_function(t_max,a0,c1a,c2a,b0,c1b,c2b,beta,alpha,delta,zetamin,zetamax,zeta_start,h,f,q,gamma,init_pop,strain_total,index_start,res0,nevol,version);
ZETA=[ZETA;ZETAnew];

% Plot the evolutionary trajectory:
subplot(2,1,2)
aa0=1;
ZETA0=log10(ZETA);
ZETA0(ZETA0<-aa0)=-aa0;
ZETA0=(ZETA0+aa0)/aa0;
imagesc(ZETA0);set(gca,'ydir','normal')
map=colormap('gray');
map=flipud(map);
colormap(map);
set(gca,'xtick',1:(res0-1)/2:res0,'xticklabel',round(Zeta(1:(res0-1)/2:res0)*10000)/10000,'fontsize',12);
set(gca,'ytick',[0,1000,2000,3000,4000]);
ylabel('Evolutionary time','interpreter','latex','fontsize',16)
xlabel('Rate of onset of resistance, $\zeta$','interpreter','latex','fontsize',16)
pbaspect([1,1.5,1])
hold on

% Plot the analytically determined value of the singular strategy:
singstrat=58/75;
singstratval=1+(singstrat/(1/(res0-1)));
singstrat_x=singstratval+zeros(4000,1);
singstrat_y=linspace(1,4000,4000);
plot(singstrat_x,singstrat_y,'--r','linewidth',2)
text(1,3900,'B','fontsize',16)
