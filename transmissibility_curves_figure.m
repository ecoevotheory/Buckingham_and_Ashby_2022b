% This code plots changes in the rate of onset of resistance (zeta) as 
% disease transmissibility (beta) varies.

%% Reproduction trade-off curves

% Define parameters:
version=3;
a0=5;
c2a=-2;
b0=0.1;
c2b=-2;
f=0.5;
h=0.75;
alpha=0;
delta=0.25;
beta=0.5;
q=1;
gamma=0;

% Set up parameters and vectors to be used later:
xmax=10; 
varvec=linspace(0,xmax,51);
ymax=0.75; 

% First curve (low costs):
c1a=0.3;
c1b=0.3;

% Set up vectors to use later:
CSS1=NaN(length(varvec),2);

% For each value of transmissibility, determine the singular strategies:
parfor i=1:length(varvec)
    
    beta=varvec(i);
    [CSSvec1,repvec1,BPvec1,GOEvec1,othervec1]=find_singstrats_function(a0,c1a,c2a,b0,c1b,c2b,f,h,beta,alpha,delta,q,gamma,version,ymax);

    % Format the output:
    CSS1temp=NaN(1,2);
    if ~isempty(CSSvec1)
        CSS1temp=[i,CSSvec1(1)];
    end
    CSS1(i,:)=CSS1temp;
    
    % We are expecting a single CSS for these particular parameter values.
    % These messages appear if this is not the case:
    if length(CSSvec1)>1
        disp("More than one CSS when i=" +i)
    end
    if ~isempty(repvec1)
        disp("Repeller when i=" +i)
    end
    if ~isempty(BPvec1)
        disp("Branching point when i=" +i)
    end
    if ~isempty(GOEvec1)
        disp("Garden of Eden when i=" +i)
    end
    if ~isempty(othervec1)
        disp("Unclassified singular strategy when i=" +i)
    end

end

% Format data for plotting:
if sum(~isnan(CSS1),'all')~=0
    xvec1=varvec(CSS1(~isnan(CSS1(:,1)),1));
    yvec1=CSS1(~isnan(CSS1(:,2)),2);
else
    xvec1=NaN(length(CSS1(:,1)),1);
    yvec1=NaN(length(CSS1(:,1)),1);
end

% Second curve (high costs):
c1a=0.8;
c1b=0.8;

CSS2=NaN(length(varvec),2);

% For each value of transmissibility, determine singular strategies:
parfor i=1:length(varvec)
    
    beta=varvec(i);
    [CSSvec2,repvec2,BPvec2,GOEvec2,othervec2]=find_singstrats_function(a0,c1a,c2a,b0,c1b,c2b,f,h,beta,alpha,delta,q,gamma,version,ymax);

    % Format the output:
    CSS2temp=NaN(1,2);
    if ~isempty(CSSvec2)
        CSS2temp=[i,CSSvec2(1)];
    end
    CSS2(i,:)=CSS2temp;
    
    % We expect a single CSS for these particular parameter values. These
    % messages appear if this is not the case:
    if length(CSSvec2)>1
        disp("More than one CSS when i=" +i)
    end
    if ~isempty(repvec2)
        disp("Repeller when i=" +i)
    end
    if ~isempty(BPvec2)
        disp("Branching point when i=" +i)
    end
    if ~isempty(GOEvec2)>1
        disp("Garden of Eden when i=" +i)
    end
    if ~isempty(othervec2)
        disp("Unclassified singular strategy when i=" +i)
    end

end

% Format output for plotting:
if sum(~isnan(CSS2),'all')~=0
    xvec2=varvec(CSS2(~isnan(CSS2(:,1)),1));
    yvec2=CSS2(~isnan(CSS2(:,2)),2);
else
    xvec2=NaN(length(CSS2(:,1)),1);
    yvec2=NaN(length(CSS2(:,1)),1);
end


%% Mortality trade-off curves

version=4;

% First curve (low costs):
c1a=0.3;
c1b=0.3;

% Set up vectors to use later:
CSS1=NaN(length(varvec),2);

% For each value of transmissibility, determine the singular strategies:
parfor i=1:length(varvec)
    
    beta=varvec(i);
    [CSSvec1,repvec1,BPvec1,GOEvec1,othervec1]=find_singstrats_function(a0,c1a,c2a,b0,c1b,c2b,f,h,beta,alpha,delta,q,gamma,version,ymax);

    % Format the output:
    CSS1temp=NaN(1,2);
    if ~isempty(CSSvec1)
        CSS1temp=[i,CSSvec1(1)];
    end
    CSS1(i,:)=CSS1temp;
    
    % We are expecting a single CSS for these particular parameter values.
    % These messages appear if this is not the case:
    if length(CSSvec1)>1
        disp("More than one CSS when i=" +i)
    end
    if ~isempty(repvec1)
        disp("Repeller when i=" +i)
    end
    if ~isempty(BPvec1)
        disp("Branching point when i=" +i)
    end
    if ~isempty(GOEvec1)
        disp("Garden of Eden when i=" +i)
    end
    if ~isempty(othervec1)
        disp("Unclassified singular strategy when i=" +i)
    end

end

% Format data for plotting:
if sum(~isnan(CSS1),'all')~=0
    xvec3=varvec(CSS1(~isnan(CSS1(:,1)),1));
    yvec3=CSS1(~isnan(CSS1(:,2)),2);
else
    xvec3=NaN(length(CSS1(:,1)),1);
    yvec3=NaN(length(CSS1(:,1)),1);
end

% Second curve (high costs):
c1a=0.8;
c1b=0.8;

CSS2=NaN(length(varvec),2);

% For each value of transmissibility, determine singular strategies:
parfor i=1:length(varvec)
    
    beta=varvec(i);
    [CSSvec2,repvec2,BPvec2,GOEvec2,othervec2]=find_singstrats_function(a0,c1a,c2a,b0,c1b,c2b,f,h,beta,alpha,delta,q,gamma,version,ymax);

    % Format the output:
    CSS2temp=NaN(1,2);
    if ~isempty(CSSvec2)
        CSS2temp=[i,CSSvec2(1)];
    end
    CSS2(i,:)=CSS2temp;
    
    % We expect a single CSS for these particular parameter values. These
    % messages appear if this is not the case:
    if length(CSSvec2)>1
        disp("More than one CSS when i=" +i)
    end
    if ~isempty(repvec2)
        disp("Repeller when i=" +i)
    end
    if ~isempty(BPvec2)
        disp("Branching point when i=" +i)
    end
    if ~isempty(GOEvec2)>1
        disp("Garden of Eden when i=" +i)
    end
    if ~isempty(othervec2)
        disp("Unclassified singular strategy when i=" +i)
    end

end

% Format output for plotting:
if sum(~isnan(CSS2),'all')~=0
    xvec4=varvec(CSS2(~isnan(CSS2(:,1)),1));
    yvec4=CSS2(~isnan(CSS2(:,2)),2);
else
    xvec4=NaN(length(CSS2(:,1)),1);
    yvec4=NaN(length(CSS2(:,1)),1);
end

%% Make the plot

% Colours for plotting:
blue=1/255*[44,123,182];
red=1/255*[215,25,28];
orange=1/255*[253,174,97];

plot(xvec1,yvec1,'-','color',blue,'markersize',10,'linewidth',3)
hold on
plot(xvec2,yvec2,':','color',blue,'markersize',10,'linewidth',3)
hold on
plot(xvec3,yvec3,'-','color',red,'markersize',10,'linewidth',3);
hold on
plot(xvec4,yvec4,':','color',red,'markersize',10,'linewidth',3);
set(gca,'xtick',[0,2,4,6,8,10],'fontsize',14)
set(gca,'ytick',[0,0.25,0.5,0.75],'fontsize',14)
ylim([0,ymax])
xlim([0,xmax])
xlabel('Disease transmissibility, $\beta$','interpreter','latex','fontsize',18)
ylabel('Rate of onset of resistance, $\zeta$','interpreter','latex','fontsize',18)
axis square