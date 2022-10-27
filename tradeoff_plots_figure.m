% This code sketches different trade-off functions and loads the model 
% schematic.

%% Reproduction trade-off

% Set up variables and parameters:
syms a(zeta)
a0=5;

% Here are some example parameter values in the case of a strong,
% decelerating trade-off:
c1a=0.5;
c2a=2;

% Define the trade-off function:
a(zeta)=a0*(1-(c1a*(1-exp(-c2a*zeta))/(1-exp(-c2a))));

% Create vectors of the values of the trade-off function:
xvec=zeros(100,1);
avec1=zeros(100,1);
for i=1:100
    xvec(i)=i/100;
    avec1(i)=a(xvec(i));
end

% Now do the same for a weak, decelerating trade-off:
c1a=0.25;
c2a=2;
a(zeta)=a0*(1-(c1a*(1-exp(-c2a*zeta))/(1-exp(-c2a))));
avec2=zeros(100,1);
for i=1:100
    avec2(i)=a(xvec(i));
end

% Now do the same for a strong, accelerating trade-off:
c1a=0.5;
c2a=-1;
a(zeta)=a0*(1-(c1a*(1-exp(-c2a*zeta))/(1-exp(-c2a))));
avec3=zeros(100,1);
for i=1:100
    avec3(i)=a(xvec(i));
end

% Now do the same for a weak, accelerating trade-off:
c1a=0.25;
c2a=-1;
a(zeta)=a0*(1-(c1a*(1-exp(-c2a*zeta))/(1-exp(-c2a))));
avec4=zeros(100,1);
for i=1:100
    avec4(i)=a(xvec(i));
end

% Plot the trade-offs:
firstcolour=1/255*[215,48,39];
secondcolour=1/255*[69,117,180];
subplot(2,3,3)
plot(xvec,avec1,'Color',firstcolour,'Linewidth',2)
hold on
plot(xvec,avec2,'Color',secondcolour,'Linewidth',2)
hold on
plot(xvec,avec3,'Color',firstcolour,'LineStyle','--','Linewidth',2) 
hold on
plot(xvec,avec4,'Color',secondcolour,'LineStyle','--','Linewidth',2)
xlabel('Rate of onset of resistance, $\zeta$','interpreter','latex','fontsize',16)
ylabel('Reproduction rate, $a(\zeta)$','interpreter','latex','fontsize',16)
ylim([2,5])
xlim([0,1])
text(0.87,4.75,'B','fontsize',24)

%%  Mortality trade-off

% Set up variables and parameters:
syms b(zeta)
b0=1;

% Here are some example parameter values in the case of a strong, 
% decelerating trade-off:
c1b=0.5;
c2b=2;

% Define the trade-off function:
b(zeta)=b0*(1+(c1b*(1-exp(-c2b*zeta))/(1-exp(-c2b))));

% Create vectors of the values of the trade-off function:
xvec=zeros(100,1);
bvec1=zeros(100,1);
for i=1:100
    xvec(i)=i/100;
    bvec1(i)=b(xvec(i));
end

% Now do the same for a weak, decelerating trade-off:
c1b=0.25;
c2b=2;
b(zeta)=b0*(1+(c1b*(1-exp(-c2b*zeta))/(1-exp(-c2b))));
bvec2=zeros(100,1);
for i=1:100
    bvec2(i)=b(xvec(i));
end

% Now do the same for a strong, accelerating trade-off:
c1b=0.5;
c2b=-1;
b(zeta)=b0*(1+(c1b*(1-exp(-c2b*zeta))/(1-exp(-c2b))));
bvec3=zeros(100,1);
for i=1:100
    bvec3(i)=b(xvec(i));
end

% Now do the same for a weak, accelerating trade-off:
c1b=0.25;
c2b=-1;
b(zeta)=b0*(1+(c1b*(1-exp(-c2b*zeta))/(1-exp(-c2b))));
bvec4=zeros(100,1);
for i=1:100
    bvec4(i)=b(xvec(i));
end

% Plot the trade-offs:
firstcolour=1/255*[215,48,39];
secondcolour=1/255*[69,117,180];
thirdcolour=1/255*[217,95,2];
fourthcolour=1/255*[27,158,119];
subplot(2,3,6)
plot(xvec,bvec1,'Color',firstcolour,'Linewidth',2)
hold on
plot(xvec,bvec2,'Color',secondcolour,'Linewidth',2)
hold on
plot(xvec,bvec3,'Color',firstcolour,'LineStyle','--','Linewidth',2) 
hold on
plot(xvec,bvec4,'Color',secondcolour,'LineStyle','--','Linewidth',2)
xlabel('Rate of onset of resistance, $\zeta$','interpreter','latex','fontsize',16)
ylabel({'Mortality rate, $b(\zeta)$'},'interpreter','latex','fontsize',16)
ylim([1,1.6])
xlim([0,1])
text(0.87,1.55,'C','fontsize',24)

%% Include model schematic image

s=subplot(2,3,[1,2,4,5]);
A=imread('Figure 1 Model Schematic.jpg');
B=imresize(A,2);
image(B)
set(gca,'xtick',[])
set(gca,'ytick',[])
s.Position=s.Position+[-0.01,0,0.01,0];
text(1250,50,'A','fontsize',24)
