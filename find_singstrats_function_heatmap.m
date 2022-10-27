function [CSSvec,other,disprev,resprop]=find_singstrats_function_heatmap(a0,c1a,c2a,b0,c1b,c2b,f,h,beta,alpha,delta,q,gamma,version,ymax)

% This function calculates singular strategies for a specific set of
% parameter values. It determines the CSS's and whether or not there are
% any other singular strategies. It also calculates the disease prevalence
% and proportion of hosts which are resistant (at the singular strategy).

% Set up parameters:
resolution=1e5;
ymin=1/resolution;

% Find and classify all singular strategies:
[CSSvec,repvec,BPvec,GOEvec,othervec]=find_singstrats_function2(a0,c1a,c2a,b0,c1b,c2b,f,h,beta,alpha,delta,q,gamma,version,ymax,ymin,resolution);

% Format output:
if isempty(CSSvec)
    CSSvec=[];
end
if isempty(repvec)
    repvec=[];
end
if isempty(BPvec)
    BPvec=[];
end
if isempty(GOEvec)
    GOEvec=[];
end
if isempty(othervec)
    othervec=[];
end

% Now, we will use simulations to classify the evolutionary stability of
% the singular strategies in the cases where we have constant costs during 
% the resistant stage (versions 1 and 2). We know analytically that E=0 in
% these cases and so evolutionary stability cannot be determined
% analytically.
if version==1 || version==2
    
    % Take all convergence stable singular strategies and determine their
    % evolutionary stability:
    reclassify1=[CSSvec;BPvec];
    CSSvectemp=NaN(length(reclassify1));
    BPvectemp=NaN(length(reclassify1));
    for i=1:length(reclassify1)
        if reclassify1(i)<0.5
            zetamin=0;
        else
            zetamin=reclassify1(i)-0.5;
        end
        zetamax=reclassify1(i)+0.5;
        if reclassify1(i)>0.2
            zetastart=reclassify1(i)-0.1;
        else
            zetastart=reclassify1(i)/2;
        end
        zetaSS=reclassify1(i);
        [isBP,~,~]=classification_simulation_function(a0,c1a,c2a,b0,c1b,c2b,f,h,beta,alpha,delta,q,gamma,zetamin,zetamax,zetastart,version,zetaSS);
        % Record whether the singular strategy is a branching point or a
        % CSS:
        if isBP==1
            BPvectemp(i)=reclassify1(i);
        else 
            CSSvectemp(i)=reclassify1(i);
        end  
    end
    CSSvectemp(isnan(CSSvectemp))=[];
    CSSvec=CSSvectemp;
    BPvectemp(isnan(BPvectemp))=[];
    BPvec=BPvectemp;
    
    % We assume that there are no 'Gardens of Eden' as they are very rare.
    % Also, to almost all intents and purposes, they are the same as
    % repellers.
    repvec=[repvec;GOEvec];
    GOEvec=[];
    
    % We take all singular strategies which were not previously classified
    % and determine their convergence and evolutionary stability through
    % simulations:
    reclassify2=othervec;
    CSSvectemp=[CSSvec;NaN(length(reclassify2))];
    BPvectemp=[BPvec;NaN(length(reclassify2))];
    repvectemp=[repvec;NaN(length(reclassify2))];
    for i=1:length(reclassify2)
        if reclassify2(i)<0.5
            zetamin=0;
        else
            zetamin=reclassify2(i)-0.5;
        end
        zetamax=reclassify2(i)+0.5;
        if reclassify2(i)>0.2
            zetastart=reclassify2(i)-0.1;
        else
            zetastart=reclassify2(i)/2;
        end
        zetaSS=reclassify2(i);
        [isBP,isCSS,~]=classification_simulation_function(a0,c1a,c2a,b0,c1b,c2b,f,h,beta,alpha,delta,q,gamma,zetamin,zetamax,zetastart,version,zetaSS);
        % Record whether the singular strategy is a branching point, CSS or
        % repeller:
        if isBP==1
            BPvectemp(length(BPvec)+i)=reclassify2(i);
        elseif isCSS==1
            CSSvectemp(length(CSSvec)+i)=reclassify2(i);
        else
            repvectemp(length(repvec)+i)=reclassify2(i);
        end
        
    end
    CSSvectemp(isnan(CSSvectemp))=[];
    CSSvec=CSSvectemp;
    BPvectemp(isnan(BPvectemp))=[];
    BPvec=BPvectemp;
    repvectemp(isnan(repvectemp))=[];
    repvec=repvectemp;
    othervec=[];
end
    
% Two singular strategies with the same stability within one 'step' of
% each other are the same singular strategy. When more than two singular
% strategies of the same type are found very close to each other, they are
% probably all one singular strategy, but its value may not be found
% accurately by averaging over them (so a warning in displayed).
if length(CSSvec)>1
    for i=2:length(CSSvec)
        for j=1:i
            if abs(CSSvec(i)-CSSvec(j))<=(ymax-ymin +1/resolution)/resolution
                CSSvec(i)=NaN;
            end
        end
    end
    if sum(isnan(CSSvec))>1
        disp("Possibly inaccurate singular strategies at version=" +version)
    end
    CSSvec(isnan(CSSvec))=[];
end
if length(repvec)>1
    for i=2:length(repvec)
        for j=1:i
            if abs(repvec(i)-repvec(j))<=(ymax-ymin +1/resolution)/resolution
                repvec(i)=NaN;
            end
        end
    end
    if sum(isnan(repvec))>1
        disp("Possibly inaccurate singular strategies at version=" +version)
    end
    repvec(isnan(repvec))=[];
end
if length(BPvec)>1
    for i=2:length(BPvec)
        for j=1:i
            if abs(BPvec(i)-BPvec(j))<=(ymax-ymin +1/resolution)/resolution
                BPvec(i)=NaN;
            end
        end
    end
    if sum(isnan(BPvec))>1
        disp("Possibly inaccurate singular strategies at version=" +version)
    end
    BPvec(isnan(BPvec))=[];
end
if length(GOEvec)>1
    for i=2:length(GOEvec)
        for j=1:i
            if abs(GOEvec(i)-GOEvec(j))<=(ymax-ymin +1/resolution)/resolution
                GOEvec(i)=NaN;
            end
        end
    end
    if sum(isnan(GOEvec))>1
        disp("Possibly inaccurate singular strategies at version=" +version)
    end
    GOEvec(isnan(GOEvec))=[];
end
if length(othervec)>1
    for i=2:length(othervec)
        for j=1:i
            if abs(othervec(i)-othervec(j))<=(ymax-ymin +1/resolution)/resolution
                othervec(i)=NaN;
            end
        end
    end
    if sum(isnan(othervec))>1
        disp("Possibly inaccurate singular strategies at version=" +version)
    end
    othervec(isnan(othervec))=[];
end

% Calculate disease prevalence and resistant proportion:
if isempty(BPvec) && isempty(GOEvec) && isempty(repvec) && isempty(othervec) && length(CSSvec)==1
    other=0;
    
    % Define trade-offs for each version of the model:
    zeta=CSSvec;
    if version==1
        bs=b0;
        as=a0;
        bRs=b0;
        aRs=a0*h;
    elseif version==2
        bs=b0;
        as=a0;
        bRs=b0*(1+delta);
        aRs=a0;
    elseif version==3
        as=a0*(1-((c1a*(1-exp(-c2a*zeta)))/(1-exp(-c2a))));
        aRs=a0*(1-((c1a*(1-exp(-c2a*zeta)))/(1-exp(-c2a))));
        bs=b0;
        bRs=b0;
    elseif version==4
        bs=b0*(1+((c1b*(1-exp(-c2b*zeta)))/(1-exp(-c2b))));
        bRs=b0*(1+((c1b*(1-exp(-c2b*zeta)))/(1-exp(-c2b))));
        as=a0;
        aRs=a0;
    elseif version==5
        as=a0*(1-((c1a*(1-exp(-c2a*zeta)))/(1-exp(-c2a))));
        aRs=a0;
        bs=b0;
        bRs=b0;
    elseif version==6
        bs=b0*(1+((c1b*(1-exp(-c2b*zeta)))/(1-exp(-c2b))));
        bRs=b0;
        as=a0;
        aRs=a0;
    end
    
    % Define useful parameters:
    A=1+alpha;
    cs=bs*A+gamma;
    S=cs/beta;
    R=(zeta/bRs)*S;
    
    % Find the density of infected individuals (I) for different values
    % of zeta:
    quad_A=as*q*f*beta*bRs;
    quad_B=-(as*(1-q*R-q*S)*f+gamma)*beta*bRs+as*q*cs*bRs+aRs*q*zeta*cs+beta*cs*bRs;
    quad_C=-as*(1-q*S-q*R)*cs*bRs-aRs*(1-q*S-q*R)*zeta*cs+(bs+zeta)*cs*bRs;
    if quad_B^2-4*quad_A*quad_C<0
        Ival=NaN;
    else
        Ival=(-quad_B+sqrt(quad_B^2-4*quad_A*quad_C))/(2*quad_A);
    end
    if quad_A==0
        Ival=-quad_C/quad_B;
    end
    if Ival<=0
        Ival=NaN;
    end
    
    % Define the total density, disease prevalence and resistant
    % proportion:
    N=S+Ival+R;
    disprev=Ival/N;
    resprop=R/N;
    
else
    % We are only considering cases where the only singular strategy is a
    % CSS. When this is not the case, this function generates the output
    % 'other=1'. 
    other=1;
    CSSvec=NaN;
    disprev=NaN;
    resprop=NaN;
end

end