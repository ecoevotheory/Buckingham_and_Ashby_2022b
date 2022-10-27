% This code plots the singular value of zeta (in the constant costs case)
% as different parameters are varied.

%% Varying f

% Define other parameters
a0=5;
b0=0.1;
q=1;
h=0.75;
delta=0;
alpha=0;
beta=1;

% Set up vectors to use later:
fvec=NaN(100,1);
zetastarvec=NaN(100,1);

% For each value of f, determine the value of the singular strategy:
for i=1:100
    f=i/100;
    fvec(i)=f;
    p=(1+delta)/h;
    zetastarvec(i)=b0*(1+delta)*(((p-1)/(f*p-1-alpha))+((beta*(a0-b0*p))/(a0*q*b0*(1+alpha)))-1);
    if h*(1+alpha)<=f*(1+delta) || zetastarvec(i)<0
        zetastarvec(i)=0;
    end
end

% Plot the curve:
subplot(2,2,1)
plot(1-fvec,zetastarvec,'linewidth',2)
xlabel('Sterility virulence, $1-f$','interpreter','latex')
ylabel('Rate of onset of resistance, $\zeta^*$','interpreter','latex')
ylim([0,1])
axis square
text(0.03,0.95,'A')
set(gca,'ytick',[0,0.5,1])

%% Varying alpha

% Define other parameters:
a0=5;
b0=0.1;
q=1;
h=1;
f=1;
delta=0.25;
beta=1;

% Set up vectors to use later:
alphavec=NaN(100,1);
zetastarvec=NaN(100,1);

% For each value of alpha, determine the value of the singular strategy:
for i=1:100
    alpha=i/100;
    alphavec(i)=alpha;
    p=(1+delta)/h;
    zetastarvec(i)=b0*(1+delta)*(((p-1)/(f*p-1-alpha))+((beta*(a0-b0*p))/(a0*q*b0*(1+alpha)))-1);
    if h*(1+alpha)<=f*(1+delta) || zetastarvec(i)<0
        zetastarvec(i)=0;
    end
end

% Plot the curve:
subplot(2,2,2)
plot(alphavec,zetastarvec,'linewidth',2)
xlabel('Mortality virulence, $\alpha$','interpreter','latex')
ylabel('Rate of onset of resistance, $\zeta^*$','interpreter','latex')
ylim([0,1])
axis square
text(0.03,0.95,'B')
set(gca,'ytick',[0,0.5,1])

%% Varying h

% Define other parameters:
a0=5;
b0=0.1;
q=1;
delta=0;
f=0.5;
alpha=0;
beta=1;

% Set up vectors to use later:
hvec=NaN(100,1);
zetastarvec=NaN(100,1);

% For each value of h, determine the value of the singular strategy:
for i=1:100
    h=i/100;
    hvec(i)=h;
    p=(1+delta)/h;
    zetastarvec(i)=b0*(1+delta)*(((p-1)/(f*p-1-alpha))+((beta*(a0-b0*p))/(a0*q*b0*(1+alpha)))-1);
    if h*(1+alpha)<=f*(1+delta) || zetastarvec(i)<0
        zetastarvec(i)=0;
    end
end

% Plot the curve:
subplot(2,2,3)
plot(1-hvec,zetastarvec,'linewidth',2)
xlabel('Fecundity cost of resistance, $1-h$','interpreter','latex')
ylabel('Rate of onset of resistance, $\zeta^*$','interpreter','latex')
ylim([0,1])
axis square
text(0.03,0.95,'C')
set(gca,'ytick',[0,0.5,1])

%% Varying delta

% Define other parameters:
a0=5;
b0=0.1;
q=1;
h=1;
f=1;
alpha=1;
beta=1;

% Set up vectors to use later:
deltavec=NaN(100,1);
zetastarvec=NaN(100,1);

% For each value of delta, determine the value of the singular strategy:
for i=1:100
    delta=i/100;
    deltavec(i)=delta;
    p=(1+delta)/h;
    zetastarvec(i)=b0*(1+delta)*(((p-1)/(f*p-1-alpha))+((beta*(a0-b0*p))/(a0*q*b0*(1+alpha)))-1);
    if h*(1+alpha)<=f*(1+delta) || zetastarvec(i)<0
        zetastarvec(i)=0;
    end
end

% Plot the curve:
subplot(2,2,4)
plot(deltavec,zetastarvec,'linewidth',2)
xlabel('Mortality cost of resistance, $\delta$','interpreter','latex')
ylabel('Rate of onset of resistance, $\zeta^*$','interpreter','latex')
ylim([0,1])
axis square
text(0.03,0.95,'D')
set(gca,'ytick',[0,0.5,1])