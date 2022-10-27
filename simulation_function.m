function [zeta_end,end_pop,strain_total,index_end,ZETA,DISPREV,RESPREV,NVEC] = simulation_function(t_max,a0,c1a,c2a,b0,c1b,c2b,beta,alpha,delta,zetamin,zetamax,zeta_start,h,f,q,gamma,init_pop,strain_total,index_start,res0,nevol,version)

% This function runs a fixed number of evolutionary timesteps for the
% evolution of zeta.

% Set up vectors and parameters:
eqtol = 1e-3;
exttol = 1e-5;
Zeta = linspace(zetamin,zetamax,res0);
ZETA = zeros(nevol,res0);
DISPREV = zeros(nevol,1);
RESPREV = zeros(nevol,1);
NVEC = zeros(nevol,1);

% Initial conditions:
zeta_current = zeta_start;
index_current = index_start;

% Each value of ievol is one evolutionary timestep.
for ievol=1:nevol
    disp(ievol)
    
    % Define the trade-offs for each version of the model. Determine the
    % ecological equilibrium of the current population.
    if version==1
        a=a0;
        b=b0;
        aR=a0*h;
        bR=b0;
        [~,S1,I1,R1,~] = Simulation_constantcost_function(t_max,a,beta,b,aR,bR,zeta_current,q,f,alpha,gamma,eqtol,init_pop,strain_total);
    elseif version==2
        a=a0;
        b=b0;
        aR=a0;
        bR=b0*(1+delta);
        [~,S1,I1,R1,~] = Simulation_constantcost_function(t_max,a,beta,b,aR,bR,zeta_current,q,f,alpha,gamma,eqtol,init_pop,strain_total);
    elseif version==3
        b=b0;
        bR=b0;
        [~,S1,I1,R1,~] = Simulation_aTradeoff_function(t_max,a0,c1a,c2a,beta,b,bR,zeta_current,q,f,alpha,gamma,eqtol,init_pop,strain_total);
    elseif version==4
        a=a0;
        aR=a0;
        [~,S1,I1,R1,~] = Simulation_bTradeoff_function(t_max,a,beta,b0,aR,c1b,c2b,zeta_current,q,f,alpha,gamma,eqtol,init_pop,strain_total);
    elseif version==5
        b=b0;
        aR=a0;
        bR=b0;
        [~,S1,I1,R1,~] = Simulation_earlyaTradeoff_function(t_max,a0,c1a,c2a,beta,b,aR,bR,zeta_current,q,f,alpha,gamma,eqtol,init_pop,strain_total);
    elseif version==6
        a=a0;
        aR=a0;
        bR=b0;
        [~,S1,I1,R1,~] = Simulation_earlybTradeoff_function(t_max,a,beta,b0,aR,bR,c1b,c2b,zeta_current,q,f,alpha,gamma,eqtol,init_pop,strain_total);
    end
        
    % Re-format this output:
    S=zeros(strain_total,1);
    I=zeros(strain_total,1);
    R=zeros(strain_total,1);
    for j=1:strain_total
        S(j,1) = S1(end,j);
        I(j,1) = I1(end,j);
        R(j,1) = R1(end,j);
    end
    N=S+I+R;
    
    % Remove extinct classes
    Ntotal=sum(N);
    extinct = (N/Ntotal)<exttol;
    strain_total = strain_total-sum(extinct);
    S(extinct) = [];
    I(extinct) = [];
    R(extinct) = [];
    N(extinct) = [];
    index_current(extinct) = [];
    zeta_current(extinct) = [];
    Ntotal=sum(N);
    
    % Proportion of individuals of each zeta strain
    ZETA(ievol,index_current) = N/Ntotal;
    % Proportion of individuals who have the disease
    DISPREV(ievol) = sum(I)/Ntotal;
    % Proportion of individuals who are resistant
    RESPREV(ievol) = sum(R)/Ntotal;
    
    % Introduce rare mutant:
    weightedprob = N/Ntotal;
    cumsum1 = cumsum(weightedprob);
    r1 = rand*cumsum1(end);
    mutator_loc = (find(r1<cumsum1,1));
    mutator = index_current(mutator_loc);
 
    if(mutator==1) % Mutate up
        mutant = mutator+1;
    elseif(mutator==res0) % Mutate down
        mutant = mutator-1;
    else
        if(rand>0.5) % Mutate up
            mutant = mutator+1;
        else % Mutate down
            mutant = mutator-1;
        end
    end
    
    if (~ismember(mutant,index_current)) % New strain
        zeta_current_new=NaN(length(zeta_current)+1,1);
        index_current_new=NaN(length(zeta_current)+1,1);
        Snew=NaN(length(zeta_current)+1,1);
        Inew=NaN(length(zeta_current)+1,1);
        Rnew=NaN(length(zeta_current)+1,1);
        for i=1:length(zeta_current)
            zeta_current_new(i)=zeta_current(i);
            index_current_new(i)=index_current(i);
            Snew(i,1)=S(i,1);
            Inew(i,1)=I(i,1);
            Rnew(i,1)=R(i,1);
        end
        strain_total=strain_total+1;
        zeta_current_new(end) = Zeta(mutant);
        index_current_new(end) = mutant;
        Snew(end,1)=S(mutator_loc,1)/10;
        Inew(end,1)=I(mutator_loc,1)/10;
        Rnew(end,1)=R(mutator_loc,1)/10;
        
        zeta_current=zeta_current_new;
        index_current=index_current_new;
        S=Snew;
        I=Inew;
        R=Rnew;
    end
    
    % Update population for the start of the next evolutionary timestep:
    N=S+I+R;
    init_pop=zeros(1,3*strain_total);
    for i=1:strain_total
        init_pop(3*i)=R(i);
        init_pop(3*i-1)=I(i);
        init_pop(3*i-2)=S(i);
    end
    Ntotal=sum(N);
    if Ntotal==0
        disp('Host has gone extinct')
        return
    end
    NVEC(ievol)=Ntotal;
    
end

% Create output:
zeta_end=zeta_current;
end_pop=init_pop;
index_end=index_current;

end