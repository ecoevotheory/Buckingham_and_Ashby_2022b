% This code simulates the evolution of the rate of onset of resistance
% (zeta) for scenarios (3) and (5) to show branching.

% Define Parameters:
a0=5;
c1a=0.1;
b0=0.1;
c1b=0.1;
f=0.5;
beta=0.5;
c2b=2;
alpha=0;
q=1;
gamma=0;
delta=0.25;
h=0.75;
t_max=100;
zetamin=0;
zetamax=5;
zetastart=0.2;
res0=101;

% First curve (decelerating trade-offs):
version=3;
c2a=5;
nevol=12000;
ZETA=[];

% Set up initial conditions:
strain_total = 1;
init_pop=[0.1,0.1,0];
Zeta = linspace(zetamin,zetamax,res0);
initial = find(Zeta>=zetastart,1);
zeta_start = Zeta(initial);
index_start = initial;

% Allow zeta to evolve:
[~,~,~,~,ZETAnew,~,~,~] = simulation_function(t_max,a0,c1a,c2a,b0,c1b,c2b,beta,alpha,delta,zetamin,zetamax,zeta_start,h,f,q,gamma,init_pop,strain_total,index_start,res0,nevol,version);
ZETA=[ZETA;ZETAnew];

% Second curve (accelerating trade-offs):
nevol=20000;
version=5;
c2a=1;
ZETA2=[];

% Allow zeta to evolve:
[~,~,~,~,ZETAnew,~,~,~] = simulation_function(t_max,a0,c1a,c2a,b0,c1b,c2b,beta,alpha,delta,zetamin,zetamax,zeta_start,h,f,q,gamma,init_pop,strain_total,index_start,res0,nevol,version);
ZETA2=[ZETA2;ZETAnew];


%% Make the Plot
aa0=1.5;
ZETA0=log10(ZETA);
ZETA0(ZETA0<-aa0)=-aa0;
ZETA0=(ZETA0+aa0)/aa0;
subplot(2,1,1)
imagesc(ZETA0);set(gca,'ydir','normal')
map=colormap('gray');
map=flipud(map);
colormap(map);
xlabel('Rate of onset of resistance, $\zeta$','interpreter','latex','fontsize',16)
ylabel('Evolutionary time','interpreter','latex','fontsize',16)
set(gca,'xtick',1:(res0-1)/10:res0,'xticklabel',round(Zeta(1:(res0-1)/10:res0)*10000)/10000);
set(gca,'ytick',[0,3000,6000,9000,12000],'yticklabel',[0,3000,6000,9000,12000]);
xlim([1,51])
ylim([0,12000])
pbaspect([1,1.5,1])
text(4.5,11500,'A','fontsize',24)

aa0=1.5;
ZETA0=-sqrt(-log10(ZETA2));
ZETA0(ZETA0<-aa0)=-aa0;
ZETA0=(ZETA0+aa0)/aa0;
subplot(2,1,2)
imagesc(ZETA0);set(gca,'ydir','normal')
map=colormap('gray');
map=flipud(map);
colormap(map);
xlabel('Rate of onset of resistance, $\zeta$','interpreter','latex','fontsize',16)
ylabel('Evolutionary time','interpreter','latex','fontsize',16)
set(gca,'xtick',1:(res0-1)/10:res0,'xticklabel',round(Zeta(1:(res0-1)/10:res0)*10000)/10000);
set(gca,'ytick',[0,5000,10000,15000,20000],'yticklabel',[0,5000,10000,15000,20000]);
xlim([1,51])
ylim([0,20000])
pbaspect([1,1.5,1])
text(4.5,19167,'B','fontsize',24)
