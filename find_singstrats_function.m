function [CSSvec,repvec,BPvec,GOEvec,othervec]=find_singstrats_function(a0,c1a,c2a,b0,c1b,c2b,f,h,beta,alpha,delta,q,gamma,version,ymax)

% This function determines the singular strategies and their stability for
% a particular set of parameter values.

% Define parameters:
resolution=1e5;
ymin=1/resolution;

% Find singular strategies:
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
% the singular strategies in the cases where costs are constant and only 
% paid during the resistant stage (versions 1 and 2). We know analytically
% that E=0 in these cases and so we can only determine evolutionary
% stability through simulations.
if version==1 || version==2
    
    % Take all of the singular strategies which are convergence stable and
    % determine their evolutionary stability:
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
        % Record whether they are branching points or CSS's:
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

    % We take all singular strategies which were previously unclassified in
    % terms of their stability and determine their convergence and
    % evolutionary stability using simulations:
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
% each other are the same singular strategy. If more than two singular
% strategies with the same stability are very close together then they are
% probably one singular strategy, but its value may not be estimated
% accurately (in which case we display a warning message).
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
    

end