function [isBP,isCSS,isother]=classification_simulation_function(a0,c1a,c2a,b0,c1b,c2b,f,h,beta,alpha,delta,q,gamma,zetamin,zetamax,zetastart,version,zetaSS)

% This function uses a simulation to determine the stability of a
% previously unclasified singular strategy.

% Set up parameters:
t_max=100;
res0=51;
nevol=10000;
ZETA=[];
DISPREV=[];
isCSS=0;
isBP=0;
isother=0;

% Initial conditions:
strain_total = 1;
init_pop=[0.1,0.1,0];
Zeta = linspace(zetamin,zetamax,res0);
initial = find(Zeta>=zetastart,1);
zeta_start = Zeta(initial);
index_start = initial;

% Allow zeta to evolve
[~,~,~,~,ZETAnew,DISPREVnew,~,~] = simulation_function(t_max,a0,c1a,c2a,b0,c1b,c2b,beta,alpha,delta,zetamin,zetamax,zeta_start,h,f,q,gamma,init_pop,strain_total,index_start,res0,nevol,version);
ZETA=[ZETA;ZETAnew];
DISPREV=[DISPREV;DISPREVnew];

% Now we classify the singular strategy based on its evolutionary
% end-point.
if DISPREV(end)==0
    % If the disease goes extinct then we cannot classify:
    isBP=0;
    isCSS=0;
    isother=1;
else
    ZETAend1=ZETA(end,:);
    ZETAend=[0 ZETAend1 0];
    [~,SSlocs]=findpeaks(ZETAend);
    singstrats=Zeta(SSlocs-1);
    % If there is more than one end-point then it must be a branching
    % point:
    if length(SSlocs)>1
        isBP=1;
    % If there is one end-point at the expected value of the singular
    % strategy then it is a CSS:
    elseif length(SSlocs)==1 && abs(singstrats-zetaSS)<3*(zetamax-zetamin)/res0
        isCSS=1;
    % Otherwise, we cannot classify the singular strategy:
    else
        isother=1;
    end
 
end

end
