% This code determines numerically whether the singular strategy for given
% sets of parameter values is convergence stable or unstable. 

% We hypothesise that it is convergence stable if and only if
% h(1+alpha)>f(1+delta). This code verifies that this hypothesis holds for
% the given sets of parameter values. 

% Set up variables:
syms zeta zetam

% Set up vectors with possible values of the different parameters:
version_vec=[1,2];
a_vec=[2,5,8];
f_vec=[0.1,0.5,0.8];
beta_vec=[0.5,1,2];
h_vec=[0.5,0.75];
b_vec=[0.1,0.4,0.7];
alpha_vec=[0.5,1,2];
delta_vec=[0.25,0.5];
q=1;

% Create vectors with all combinations of these parameter values:
all_combinations=combvec(version_vec,a_vec,f_vec,beta_vec,h_vec,b_vec,alpha_vec,delta_vec);
versionvec=all_combinations(1,:);
avec=all_combinations(2,:);
fvec=all_combinations(3,:);
betavec=all_combinations(4,:);
hvec=all_combinations(5,:);
bvec=all_combinations(6,:);
alphavec=all_combinations(7,:);
deltavec=all_combinations(8,:);

% Set up vectors to be used later:
Mvec=NaN(length(avec),1);
negvec=NaN(length(avec),1);

% For each combination of parameters, we find the singular strategy and
% determine its convergence stability:
for k=1:length(avec)
    
    % Define a particular combination of parameter values:
    version=versionvec(k);
    a=avec(k);
    if version==1
        f=fvec(k);
        alpha=0;
    end
    beta=betavec(k);
    h=hvec(k);
    b=bvec(k);
    if version==2
        alpha=alphavec(k);
        f=1;
    end
    delta=deltavec(k);
    
    % Find the endemic equilibrium:
    Sstar=b*(1+alpha)/beta;
    Rstar=zeta*Sstar/(b*(1+delta));
    A=a*f;
    B=a*f*Sstar+a*f*Rstar-a*f+beta*Sstar+a*Sstar+a*h*Rstar;
    C=b*Sstar+zeta*Sstar-a*Sstar-a*h*Rstar+a*Sstar^2+a*h*Sstar*Rstar+a*Sstar*Rstar+a*h*Rstar^2;
    Istar=(-B+sqrt(B^2-4*A*C))/(2*A);
    Nstar=Sstar+Istar+Rstar;
    
    % Define the invasion fitness:
    w=((a*(1-q*Nstar)*b*(1+alpha)*b*(1+delta)+a*(1-q*Nstar)*f*beta*Istar*b*(1+delta)+a*h*(1-q*Nstar)*b*(1+alpha)*zetam)/((beta*Istar+b+zetam)*b*(1+alpha)*b*(1+delta)))-1;
    
    % Find the fitness gradient and the second cross-derivative of the 
    % invasion fitness (which determines the convergence stability since
    % E=0):
    fitgrad=diff(w,zetam);
    Mfunction=diff(fitgrad,zeta);
    
    % Calculate the singular strategy:
    fitgrad_1var=subs(fitgrad,zetam,zeta);
    singstrat=solve(fitgrad_1var,zeta);
    if singstrat<0
        singstrat=0;
        negvec(k)=0;
    end
    
    % If no singular strategy is found, determine whether the fitness
    % gradient is always positive or always negative:
    if isempty(singstrat)
        zetaval=0.4;
        fitgradplotter=subs(fitgrad_1var,zeta,zetaval);
        if fitgradplotter>0
            negvec(k)=1;
        elseif fitgradplotter<0
            singstrat=0;
            negvec(k)=-1;
        end
    end
    
    % If there is a singular strategy, find the value of the second
    % cross-derivative at the singular strategy:
    if ~isempty(singstrat) && isnan(negvec(k))
        Mvec(k)=double(subs(Mfunction,[zeta,zetam],[singstrat,singstrat]));
    end
    
    % We expect the second cross-derivative and our original inequality to
    % have different signs (this is our hypothesis). If this fails, then a
    % warning message is displayed:
    if Mvec(k)>0 && h*(1+alpha)-f*(1+delta)>0
        disp("M and inequality both positive when k=" +k)
    elseif Mvec(k)<0 && h*(1+alpha)-f*(1+delta)<0
        disp("M and inequality both negative when k=" +k)
    end

end