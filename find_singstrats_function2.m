function [CSSvec,repvec,BPvec,GOEvec,othervec]=find_singstrats_function2(a0,c1a,c2a,b0,c1b,c2b,f,h,beta,alpha,delta,q,gamma,version,ymax,ymin,resolution)

% This function finds singular strategies and determines their evolutionary
% and convergence stability.

% Set up parameters and vectors to use later:
max_singstrat_value=ymax;
min_singstrat_value=ymin;
if version~=1 && version~=2 && version~=3 && version~=4 && version~=5 && version~=6
    disp('Please select which version of the model you are using')
end
fitgradmat=NaN(resolution,3);

% First, we determine the fitness gradient for different values of zeta 
% (resident) and zetam (mutant). Note that we only care about sign changes 
% when zeta=zetam, and so we only need to find the fitness gradient when 
% zeta is close to zetam.
zetavec=linspace(min_singstrat_value,max_singstrat_value,resolution);
for i=1:resolution
    
    zeta=zetavec(i);
    
    % We will look at the fitness gradient when zetam is within one step of
    % zeta:
    if i==1
        startval=1;
        endval=2;
    elseif i==resolution
        startval=resolution-1;
        endval=resolution;
    else
        startval=i-1;
        endval=i+1;
    end
    
    for j=startval:endval
        zetam=zetavec(j);
    
        % For each version, set up the costs. The suffix 's' represents the
        % functions of resident (non-mutant) zeta.
        if version==1
            b=b0;
            a=a0;
            bR=b0;
            aR=a0*h;
            dadz=0;
            daRdz=0;
            dbdz=0;
            dbRdz=0;
            bs=b0;
            as=a0;
            bRs=b0;
            aRs=a0*h;
        elseif version==2
            b=b0;
            a=a0;
            bR=b0*(1+delta);
            aR=a0;
            dadz=0;
            daRdz=0;
            dbdz=0;
            dbRdz=0;
            bs=b0;
            as=a0;
            bRs=b0*(1+delta);
            aRs=a0;
        elseif version==3
            a=a0*(1-((c1a*(1-exp(-c2a*zetam)))/(1-exp(-c2a))));
            dadz=-a0*c1a*c2a*exp(-c2a*zetam)/(1-exp(-c2a));
            aR=a0*(1-((c1a*(1-exp(-c2a*zetam)))/(1-exp(-c2a))));
            daRdz=-a0*c1a*c2a*exp(-c2a*zetam)/(1-exp(-c2a));
            as=a0*(1-((c1a*(1-exp(-c2a*zeta)))/(1-exp(-c2a))));
            aRs=a0*(1-((c1a*(1-exp(-c2a*zeta)))/(1-exp(-c2a))));
            b=b0;
            bR=b0;
            dbdz=0;
            dbRdz=0;
            bs=b0;
            bRs=b0;
        elseif version==4
            b=b0*(1+((c1b*(1-exp(-c2b*zetam)))/(1-exp(-c2b))));
            dbdz=b0*c1b*c2b*exp(-c2b*zetam)/(1-exp(-c2b));
            bR=b0*(1+((c1b*(1-exp(-c2b*zetam)))/(1-exp(-c2b))));
            dbRdz=b0*c1b*c2b*exp(-c2b*zetam)/(1-exp(-c2b));
            bs=b0*(1+((c1b*(1-exp(-c2b*zeta)))/(1-exp(-c2b))));
            bRs=b0*(1+((c1b*(1-exp(-c2b*zeta)))/(1-exp(-c2b))));
            a=a0;
            aR=a0;
            dadz=0;
            daRdz=0;
            as=a0;
            aRs=a0;
        elseif version==5
            a=a0*(1-((c1a*(1-exp(-c2a*zetam)))/(1-exp(-c2a))));
            dadz=-a0*c1a*c2a*exp(-c2a*zetam)/(1-exp(-c2a));
            as=a0*(1-((c1a*(1-exp(-c2a*zeta)))/(1-exp(-c2a))));
            aR=a0;
            daRdz=0;
            b=b0;
            bR=b0;
            dbdz=0;
            dbRdz=0;
            aRs=a0;
            bs=b0;
            bRs=b0;
        elseif version==6
            b=b0*(1+((c1b*(1-exp(-c2b*zetam)))/(1-exp(-c2b))));
            dbdz=b0*c1b*c2b*exp(-c2b*zetam)/(1-exp(-c2b));
            bs=b0*(1+((c1b*(1-exp(-c2b*zeta)))/(1-exp(-c2b))));
            bR=b0;
            dbRdz=0;
            a=a0;
            aR=a0;
            dadz=0;
            daRdz=0;
            bRs=b0;
            as=a0;
            aRs=a0;
        end
    
        % Define useful parameters for notational convenience:
        A=1+alpha;
        c=b*A+gamma;
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
        
        if isempty(Ival) || isnan(Ival)
            S=((bRs)/(q*(bRs+zeta)))- (((bs+zeta)*bRs^2)/(q*(bRs+zeta)*(as*bRs+aRs*zeta)));
            R=zeta*S/bRs;
            Ival=0;
        end
        N=S+Ival+R;
        
        % Determine the fitness gradient from the value of I*. These
        % expressions are calculated analytically.
        term1=(aR*(1-q*N)*c)/((beta*Ival+b+zetam)*c*bR);
        term2=(a*(1-q*N)*c*bR+(a*(1-q*N)*f+gamma)*beta*Ival*bR+aR*(1-q*N)*c*zetam)/((beta*Ival+b+zetam)^2*c*bR);
        term3=((1-q*N)*c*bR+(1-q*N)*f*beta*Ival*bR)/((beta*Ival+b+zetam)*c*bR);
        term4=((1-q*N)*c*zetam)/((beta*Ival+b+zetam)*c*bR);
        t5f1=(a*(1-q*N)*A*bR+aR*(1-q*N)*A*zetam)/((beta*Ival+b+zetam)*c*bR);
        t5f2=(a*(1-q*N)*c*bR+(a*(1-q*N)*f+gamma)*beta*Ival*bR+aR*(1-q*N)*c*zetam)/((beta*Ival+b+zetam)^2*c*bR);
        t5f3=(A*a*(1-q*N)*c*bR+A*(a*(1-q*N)*f+gamma)*beta*Ival*bR+A*aR*(1-q*N)*c*zetam)/((beta*Ival+b+zetam)*c^2*bR);
        t6f1=(a*(1-q*N)*c+(a*(1-q*N)*f+gamma)*beta*Ival)/((beta*Ival+b+zetam)*c*bR);
        t6f2=(a*(1-q*N)*c*bR+(a*(1-q*N)*f+gamma)*beta*Ival*bR+aR*(1-q*N)*c*zetam)/((beta*Ival+b+zetam)*c*bR^2);
        fitgrad=(term1)-(term2)+dadz*(term3)+daRdz*(term4)+dbdz*((t5f1)-(t5f2)-(t5f3))+dbRdz*((t6f1)-(t6f2));
        
        % Record the fitness gradient:
        if j==i-1
            fitgradmat(i,1)=fitgrad;
        elseif j==i
            fitgradmat(i,2)=fitgrad;
        elseif j==i+1
            fitgradmat(i,3)=fitgrad;
        end
   
    end
end

% Now we wil use the fitness gradients to find the singular strategies and
% their stability.
% In a matrix of fitness gradients for all values of zeta and zetam, the 
% singular strategy is determine by where the main diagonal changes
% sign. The evolutionary stability is determined by the direction of the
% sign change along the row. The convergence stability is determined by the
% direction of the sign change down the column.
repvec=NaN(resolution-1,1);
BPvec=NaN(resolution-1,1);
CSSvec=NaN(resolution-1,1);
GOEvec=NaN(resolution-1,1);
othervec=NaN(resolution-1,1);
for i=2:resolution-1
    % The singular strategy could be a repeller:
    if fitgradmat(i,1)<fitgradmat(i,2) && fitgradmat(i,3)>fitgradmat(i,2) && fitgradmat(i,3)+fitgradmat(i+1,1)>fitgradmat(i,1)+fitgradmat(i-1,3) && ( (fitgradmat(i,2)>=0 && fitgradmat(i-1,2)<=0 ) || (fitgradmat(i,2)<=0 && fitgradmat(i-1,2)>=0) )
        if isnan(repvec(i-1))
            repvec(i)=zetavec(i);
        end
    % It could be a branching point:
    elseif fitgradmat(i,1)<fitgradmat(i,2) && fitgradmat(i,3)>fitgradmat(i,2) && fitgradmat(i,3)+fitgradmat(i+1,1)<fitgradmat(i,1)+fitgradmat(i-1,3) && ( (fitgradmat(i,2)>=0 && fitgradmat(i-1,2)<=0 ) || (fitgradmat(i,2)<=0 && fitgradmat(i-1,2)>=0) )
        if isnan(BPvec(i-1))
            BPvec(i)=zetavec(i);
        end
    % It could be a CSS:
    elseif fitgradmat(i,1)>fitgradmat(i,2) && fitgradmat(i,3)<fitgradmat(i,2) && fitgradmat(i,3)+fitgradmat(i+1,1)<fitgradmat(i,1)+fitgradmat(i-1,3) && ( (fitgradmat(i,2)>=0 && fitgradmat(i-1,2)<=0 ) || (fitgradmat(i,2)<=0 && fitgradmat(i-1,2)>=0) )
        if isnan(CSSvec(i-1))
            CSSvec(i)=zetavec(i);
        end
    % Or it could be a 'Garden of Eden':
    elseif fitgradmat(i,1)>fitgradmat(i,2) && fitgradmat(i,3)<fitgradmat(i,2) && fitgradmat(i,3)+fitgradmat(i+1,1)>fitgradmat(i,1)+fitgradmat(i-1,3) && ( (fitgradmat(i,2)>=0 && fitgradmat(i-1,2)<=0 ) || (fitgradmat(i,2)<=0 && fitgradmat(i-1,2)>=0) )
        if isnan(GOEvec(i-1))
            GOEvec(i)=zetavec(i);
        end
    % If the evolutionary stability condition, E, is zero, then we cannot
    % categorise the evolutionary stability:
    elseif (fitgradmat(i,2)>=0 && fitgradmat(i-1,2)<=0 ) || (fitgradmat(i,2)<=0 && fitgradmat(i-1,2)>=0)
        if isnan(othervec(i-1))
            othervec(i)=zetavec(i);
        end
    end
end
repvec(isnan(repvec))=[];
BPvec(isnan(BPvec))=[];
CSSvec(isnan(CSSvec))=[];
GOEvec(isnan(GOEvec))=[];
othervec(isnan(othervec))=[];

% If the fitness gradient is always negative then the trait will evolve to
% zero. If the fitness gradient is always positive then the trait will
% evolve to increase indefinitely.
if sum(isnan(fitgradmat),'all')==resolution*3
    % If the disease is never viable then there is no singular strategy:
    disp("No singular strategy - disease not viable")
elseif sum(fitgradmat<0,'all')+sum(isnan(fitgradmat),'all')==resolution*3
    % Fitness gradient always negative:
    CSSvec=0;
elseif sum(fitgradmat>0,'all')+sum(isnan(fitgradmat),'all')==resolution*3
    % Fitness gradient always positive:
    repvec=0;
end

end