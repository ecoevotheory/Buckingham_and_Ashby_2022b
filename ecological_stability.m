% This code shows that the endemic equilibrium is always linearly stable 
% by carrying out linear stability analysis numerically for a range of 
% parameter values.

% Set up useful quantites and parameters which are not varied:
initvec=[0.1,0.1,0.1,0.1];
t_max=1000;
maxrange=0.001;
b0=1;
q=1;

% We consider a variety of parameter values of the orders of magnitude
% considered in the paper:
version_vec=[1,2,3,4,5,6];
a0_vec=[1,3,5];
b0_vec=[0.1,0.2,0.5];
c1a_vec=0.1;
c2a_vec=[-2,2];
c1b_vec=0.1;
c2b_vec=[-2,2];
f_vec=[0.25,0.5,0.75,1];
beta_vec=[0.5,1];
alpha_vec=[0,0.5,1,2];
h_vec=[0.25,0.5,0.75];
delta_vec=[0.5,1,2];
zeta_vec=[0.1,0.3,0.5,0.7,0.9];

% Create vectors of all possible combinations of these parameter values:
all_combinations=combvec(version_vec,a0_vec,b0_vec,c1a_vec,c2a_vec,c1b_vec,c2b_vec,f_vec,beta_vec,alpha_vec,h_vec,delta_vec,zeta_vec);
versionvec=all_combinations(1,:);
a0vec=all_combinations(2,:);
b0vec=all_combinations(3,:);
c1avec=all_combinations(4,:);
c2avec=all_combinations(5,:);
c1bvec=all_combinations(6,:);
c2bvec=all_combinations(7,:);
fvec=all_combinations(8,:);
betavec=all_combinations(9,:);
alphavec=all_combinations(10,:);
hvec=all_combinations(11,:);
deltavec=all_combinations(12,:);
zetavec=all_combinations(13,:);

% For each set of parameters, we determine the endemic equilibrium
% analytically, find the Jacobian matrix of the linearised system and
% calculate its eigenvalues (which tell us linear stability).
isproblem=0;
for k=1:length(a0vec)
    disp(k)
    
    % Define parameter values:
    version=versionvec(k);
    a0=a0vec(k);
    b0=b0vec(k);
    c1a=c1avec(k);
    c2a=c2avec(k);
    c1b=c1bvec(k);
    c2b=c2bvec(k);
    f=fvec(k);
    beta=betavec(k);
    alpha=alphavec(k);
    h=hvec(k);
    delta=deltavec(k);
    zeta=zetavec(k);
    if f==1 && alpha==0
        f=0.75;
    end

    % Define trade-offs:
    if version==1
        a=a0;
        b=b0;
        aR=a0*h;
        bR=b0;
    elseif version==2
        a=a0;
        b=b0;
        aR=a0;
        bR=b0*(1+delta);
    elseif version==3
        a=a0*(1-(c1a*(1-exp(-c2a*zeta)))/(1-exp(-c2a)));
        b=b0;
        aR=a0*(1-(c1a*(1-exp(-c2a*zeta)))/(1-exp(-c2a)));
        bR=b0;
    elseif version==4
        a=a0;
        b=b0*(1+(c1b*(1-exp(-c2b*zeta)))/(1-exp(-c2b)));
        aR=a0;
        bR=b0*(1+(c1b*(1-exp(-c2b*zeta)))/(1-exp(-c2b)));
    elseif version==5
        a=a0*(1-(c1a*(1-exp(-c2a*zeta)))/(1-exp(-c2a)));
        b=b0;
        aR=a0;
        bR=b0;
    elseif version==6
        a=a0;
        b=b0*(1+(c1b*(1-exp(-c2b*zeta)))/(1-exp(-c2b)));
        aR=a0;
        bR=b0;
    end
    
    % Calculate useful parameters:
    Sstar=b*(1+alpha)/beta;
    Rstar=zeta*Sstar/bR;
    A=a*f;
    B=a*f*Sstar+a*f*Rstar-a*f+beta*Sstar+a*Sstar+aR*Rstar;
    C=b*Sstar+zeta*Sstar-a*Sstar-aR*Rstar+a*Sstar^2+aR*Sstar*Rstar+a*Sstar*Rstar+aR*Rstar^2;
    
    if C<0 % if the pathogen is viable
        
        % Determine the endemic equilibrium:
        Istar=(-B+sqrt(B^2-4*A*C))/(2*A);
        Nstar=Sstar+Istar+Rstar;
        
        % Determine the Jacobian matrix of the linearised system:
        M11=a*(1-q*Nstar)-a*q*Sstar-a*q*f*Istar-aR*q*Rstar-beta*Istar-b-zeta;
        M12=a*f-a*q*Sstar-a*q*f*Nstar-a*q*f*Istar-aR*q*Rstar-beta*Sstar;
        M13=-a*q*Sstar-a*q*f*Istar+aR-aR*q*Nstar-aR*q*Rstar;
        M21=beta*Istar;
        M22=beta*Sstar-b*(1+alpha);
        M23=0;
        M31=zeta;
        M32=0;
        M33=-bR;
        M=[M11,M12,M13;M21,M22,M23;M31,M32,M33];
        
        % Find its eigenvalues:
        eigenvalues=double(eig(M));
        if sum(real(eigenvalues)<0)~=length(eigenvalues)
            % If the eigenvalues have positive real part then the
            % equilibrium is unstable:
            disp("Eigenvalues do NOT all have negative real part at k=" +k)
            isproblem=1;
        end
    
    end
    
end

% Display this message if all of the equilibria are stable:
if isproblem==0
	disp('All equilibria are linearly stable')
end
