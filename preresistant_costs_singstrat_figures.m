% This code plots the effect of different parameters on the number and
% stability of singular strategies (in the case where costs are paid only
% before the onset of resistance). 

% Colours for plotting:
green=1/255*[27,158,119];
orange=1/255*[217,95,2];
purple=1/255*[117,112,179];
paleblue=[1,0.6,1];
midblue=[0,0.6,0.1];

%% Whole-life trade-offs, effect of c2a

% Define parameters:
a0=5;
c1a=0.1;
b0=0.1;
c1b=0.1;
c2b=2;
f=1;
h=0.75;
beta=0.5;
alpha=1;
delta=0.25;
q=1;
gamma=0;
ymax=1;
xaxislabel='$c_2^a$';
version=5;
varvec=linspace(0.1,10,100);

% Set up vectors to use later:
CSS=NaN(length(varvec),2);
rep=NaN(length(varvec),2);
BP=NaN(length(varvec),2);
GOE=NaN(length(varvec),2);
other=NaN(length(varvec),2);
disp("working on plot 1")

% For each value of c2a, we determine the singular strategies and their
% stability:
parfor i=1:length(varvec)
    
    % Set up vectors:
    CSStemp=NaN(1,2);
    reptemp=NaN(1,2);
    BPtemp=NaN(1,2);
    GOEtemp=NaN(1,2);
    othertemp=NaN(1,2);
    
    % For different combinations of trade-offs, we find the singular
    % strategies and their stability:
    c2a=varvec(i);
    [CSSvec,repvec,BPvec,GOEvec,othervec]=find_singstrats_function(a0,c1a,c2a,b0,c1b,c2b,f,h,beta,alpha,delta,q,gamma,version,ymax);

    % Format output:
    if ~isempty(CSSvec)
        CSStemp=[i,CSSvec(1)];
    end
    if ~isempty(repvec)
        reptemp=[i,repvec(1)];
    end
    if ~isempty(BPvec)
        BPtemp=[i,BPvec(1)];
    end
    if ~isempty(GOEvec)
        GOEtemp=[i,GOEvec(1)];
    end
    if ~isempty(othervec)
        othertemp=[i,othervec(1)];
    end
    
    % It is very rare to see more than one singular strategy of the same
    % type. In these cases, the following messages appear:
    if length(CSSvec)>1
        disp("More than one CSS when i=" +i)
    end
    if length(repvec)>1
        disp("More than one repeller when i=" +i)
    end
    if length(BPvec)>1
        disp("More than one BP when i=" +i)
    end
    if length(GOEvec)>1
        disp("More than one GOE when i=" +i)
    end
    if length(othervec)>1
        disp("More than one 'other' when i=" +i)
    end
    
    % Create output:
    CSS(i,:)=CSStemp;
    rep(i,:)=reptemp;
    BP(i,:)=BPtemp;
    GOE(i,:)=GOEtemp;
    other(i,:)=othertemp;

end

% Format data:
for i=1:length(varvec)
    if min([CSS(i,2),BP(i,2),rep(i,2),other(i,2),GOE(i,2)])==min(CSS(i,2),BP(i,2)) && isnan(rep(i,2))
        rep(i,2)=0;
        rep(i,1)=i;
    end
end
if sum(~isnan(CSS),'all')~=0
    xvec=varvec(CSS(~isnan(CSS(:,1)),1));
    yvec=CSS(~isnan(CSS(:,2)),2);
else
    xvec=NaN(length(CSS(:,1)),1);
    yvec=NaN(length(CSS(:,1)),1);
end
if sum(~isnan(rep),'all')~=0
    xrepvec=varvec(rep(~isnan(rep(:,1)),1));
    yrepvec=rep(~isnan(rep(:,2)),2);
else
    xrepvec=NaN(length(rep(:,1)),1);
    yrepvec=NaN(length(rep(:,1)),1);
end
if sum(~isnan(BP),'all')~=0
    xBPvec=varvec(BP(~isnan(BP(:,1)),1));
    yBPvec=BP(~isnan(BP(:,2)),2);
else
    xBPvec=NaN(length(BP(:,1)),1);
    yBPvec=NaN(length(BP(:,1)),1);
end
if sum(~isnan(GOE),'all')~=0
    xGOEvec=varvec(GOE(~isnan(GOE(:,1)),1));
    yGOEvec=GOE(~isnan(GOE(:,2)),2);
else
    xGOEvec=NaN(length(GOE(:,1)),1);
    yGOEvec=NaN(length(GOE(:,1)),1);
end
if sum(~isnan(other),'all')~=0
    xothervec=varvec(other(~isnan(other(:,1)),1));
    yothervec=other(~isnan(other(:,2)),2);
else
    xothervec=NaN(length(other(:,1)),1);
    yothervec=NaN(length(other(:,1)),1);
end
xrepvec1=NaN(length(yrepvec),1);
xrepvec2=NaN(length(yrepvec),1);
yrepvec1=NaN(length(yrepvec),1);
yrepvec2=NaN(length(yrepvec),1);
for i=1:length(yrepvec)
    if yrepvec(i)>0.1
        yrepvec1(i)=yrepvec(i);
        xrepvec1(i)=xrepvec(i);
    else
        yrepvec2(i)=yrepvec(i);
        xrepvec2(i)=xrepvec(i);
    end
end
xrepvec1(isnan(xrepvec1))=[];
xrepvec2(isnan(xrepvec2))=[];
yrepvec1(isnan(yrepvec1))=[];
yrepvec2(isnan(yrepvec2))=[];

% Create the plot:
subplot(3,3,1)
plot(xGOEvec,yGOEvec,'x','color',paleblue,'markersize',8,'linewidth',2)
axis square
hold on
plot(xBPvec,yBPvec,'--','color',purple,'markersize',8,'linewidth',3)
hold on
plot(xrepvec1,yrepvec1,':','color',orange,'markersize',8,'linewidth',3)
hold on
plot(xrepvec2,yrepvec2,':','color',orange,'markersize',8,'linewidth',3)
hold on
plot(xvec,yvec,'-','color',green,'markersize',8,'linewidth',3)
hold on
plot(xothervec,yothervec,'x','color',midblue,'markersize',8,'linewidth',2)
set(gca,'xtick',[0,10],'fontsize',14)
set(gca,'ytick',[0,1],'fontsize',14)
ylim([0,ymax])
text(0.5,0.9,'A','fontsize',16)
xlabel(xaxislabel,'interpreter','latex');
axis square

%% Whole-life trade-offs, effect of c2b

% Define parameters:
a0=5;
c1a=0.1;
c2a=2;
b0=0.1;
c1b=0.1;
f=1;
h=0.75;
beta=0.5;
alpha=1;
delta=0.25;
q=1;
gamma=0;
ymax=1;
xaxislabel='$c_2^b$';
version=6;
varvec=linspace(0.1,10,100);

% Set up vectors to use later:
CSS=NaN(length(varvec),2);
rep=NaN(length(varvec),2);
BP=NaN(length(varvec),2);
GOE=NaN(length(varvec),2);
other=NaN(length(varvec),2);
disp("working on plot 2")

% For each value of c2b, we determine the singular strategies and their
% stability:
parfor i=1:length(varvec)
    
    % Set up vectors:
    CSStemp=NaN(1,2);
    reptemp=NaN(1,2);
    BPtemp=NaN(1,2);
    GOEtemp=NaN(1,2);
    othertemp=NaN(1,2);
    
    % For different combinations of trade-offs, we find the singular
    % strategies and their stability:
    c2b=varvec(i);
    [CSSvec,repvec,BPvec,GOEvec,othervec]=find_singstrats_function(a0,c1a,c2a,b0,c1b,c2b,f,h,beta,alpha,delta,q,gamma,version,ymax);

    % Format output:
    if ~isempty(CSSvec)
        CSStemp=[i,CSSvec(1)];
    end
    if ~isempty(repvec)
        reptemp=[i,repvec(1)];
    end
    if ~isempty(BPvec)
        BPtemp=[i,BPvec(1)];
    end
    if ~isempty(GOEvec)
        GOEtemp=[i,GOEvec(1)];
    end
    if ~isempty(othervec)
        othertemp=[i,othervec(1)];
    end
    
    % It is very rare to see more than one singular strategy of the same
    % type. In these cases, the following messages appear:
    if length(CSSvec)>1
        disp("More than one CSS when i=" +i)
    end
    if length(repvec)>1
        disp("More than one repeller when i=" +i)
    end
    if length(BPvec)>1
        disp("More than one BP when i=" +i)
    end
    if length(GOEvec)>1
        disp("More than one GOE when i=" +i)
    end
    if length(othervec)>1
        disp("More than one 'other' when i=" +i)
    end
    
    % Create output:
    CSS(i,:)=CSStemp;
    rep(i,:)=reptemp;
    BP(i,:)=BPtemp;
    GOE(i,:)=GOEtemp;
    other(i,:)=othertemp;

end

% Format data:
for i=1:length(varvec)
    if min([CSS(i,2),BP(i,2),rep(i,2),other(i,2),GOE(i,2)])==min(CSS(i,2),BP(i,2)) && isnan(rep(i,2))
        rep(i,2)=0;
        rep(i,1)=i;
    end
end
if sum(~isnan(CSS),'all')~=0
    xvec=varvec(CSS(~isnan(CSS(:,1)),1));
    yvec=CSS(~isnan(CSS(:,2)),2);
else
    xvec=NaN(length(CSS(:,1)),1);
    yvec=NaN(length(CSS(:,1)),1);
end
if sum(~isnan(rep),'all')~=0
    xrepvec=varvec(rep(~isnan(rep(:,1)),1));
    yrepvec=rep(~isnan(rep(:,2)),2);
else
    xrepvec=NaN(length(rep(:,1)),1);
    yrepvec=NaN(length(rep(:,1)),1);
end
if sum(~isnan(BP),'all')~=0
    xBPvec=varvec(BP(~isnan(BP(:,1)),1));
    yBPvec=BP(~isnan(BP(:,2)),2);
else
    xBPvec=NaN(length(BP(:,1)),1);
    yBPvec=NaN(length(BP(:,1)),1);
end
if sum(~isnan(GOE),'all')~=0
    xGOEvec=varvec(GOE(~isnan(GOE(:,1)),1));
    yGOEvec=GOE(~isnan(GOE(:,2)),2);
else
    xGOEvec=NaN(length(GOE(:,1)),1);
    yGOEvec=NaN(length(GOE(:,1)),1);
end
if sum(~isnan(other),'all')~=0
    xothervec=varvec(other(~isnan(other(:,1)),1));
    yothervec=other(~isnan(other(:,2)),2);
else
    xothervec=NaN(length(other(:,1)),1);
    yothervec=NaN(length(other(:,1)),1);
end
xrepvec1=NaN(length(yrepvec),1);
xrepvec2=NaN(length(yrepvec),1);
yrepvec1=NaN(length(yrepvec),1);
yrepvec2=NaN(length(yrepvec),1);
for i=1:length(yrepvec)
    if yrepvec(i)>0.1
        yrepvec1(i)=yrepvec(i);
        xrepvec1(i)=xrepvec(i);
    else
        yrepvec2(i)=yrepvec(i);
        xrepvec2(i)=xrepvec(i);
    end
end
xrepvec1(isnan(xrepvec1))=[];
xrepvec2(isnan(xrepvec2))=[];
yrepvec1(isnan(yrepvec1))=[];
yrepvec2(isnan(yrepvec2))=[];

% Create the plot:
subplot(3,3,2)
plot(xGOEvec,yGOEvec,'x','color',paleblue,'markersize',8,'linewidth',2)
axis square
hold on
plot(xBPvec,yBPvec,'--','color',purple,'markersize',8,'linewidth',3)
hold on
plot(xrepvec1,yrepvec1,':','color',orange,'markersize',8,'linewidth',3)
hold on
plot(xrepvec2,yrepvec2,':','color',orange,'markersize',8,'linewidth',3)
hold on
plot(xvec,yvec,'-','color',green,'markersize',8,'linewidth',3)
hold on
plot(xothervec,yothervec,'x','color',midblue,'markersize',8,'linewidth',2)
set(gca,'xtick',[0,10],'fontsize',14)
set(gca,'ytick',[0,1],'fontsize',14)
ylim([0,ymax])
text(0.5,0.9,'B','fontsize',16)
xlabel(xaxislabel,'interpreter','latex');
axis square

%% Whole-life trade-offs, effect of a0

% Define parameters:
c1a=0.1;
c2a=2;
b0=0.1;
c1b=0.1;
c2b=2;
f=0.5;
h=0.75;
beta=0.5;
alpha=0;
delta=0.25;
q=1;
gamma=0;
ymax=1;
xaxislabel='$a_0$';
version=6;
varvec=linspace(0.01,1,100);

% Set up vectors to use later:
CSS=NaN(length(varvec),2);
rep=NaN(length(varvec),2);
BP=NaN(length(varvec),2);
GOE=NaN(length(varvec),2);
other=NaN(length(varvec),2);
disp("working on plot 3")

% For each value of a0, we determine the singular strategies and their
% stability:
parfor i=1:length(varvec)
    
    % Set up vectors:
    CSStemp=NaN(1,2);
    reptemp=NaN(1,2);
    BPtemp=NaN(1,2);
    GOEtemp=NaN(1,2);
    othertemp=NaN(1,2);
    
    % For different combinations of trade-offs, we find the singular
    % strategies and their stability:
    a0=varvec(i);
    [CSSvec,repvec,BPvec,GOEvec,othervec]=find_singstrats_function(a0,c1a,c2a,b0,c1b,c2b,f,h,beta,alpha,delta,q,gamma,version,ymax);

    % Format output:
    if ~isempty(CSSvec)
        CSStemp=[i,CSSvec(1)];
    end
    if ~isempty(repvec)
        reptemp=[i,repvec(1)];
    end
    if ~isempty(BPvec)
        BPtemp=[i,BPvec(1)];
    end
    if ~isempty(GOEvec)
        GOEtemp=[i,GOEvec(1)];
    end
    if ~isempty(othervec)
        othertemp=[i,othervec(1)];
    end
    
    % It is very rare to see more than one singular strategy of the same
    % type. In these cases, the following messages appear:
    if length(CSSvec)>1
        disp("More than one CSS when i=" +i)
    end
    if length(repvec)>1
        disp("More than one repeller when i=" +i)
    end
    if length(BPvec)>1
        disp("More than one BP when i=" +i)
    end
    if length(GOEvec)>1
        disp("More than one GOE when i=" +i)
    end
    if length(othervec)>1
        disp("More than one 'other' when i=" +i)
    end
    
    % Create output:
    CSS(i,:)=CSStemp;
    rep(i,:)=reptemp;
    BP(i,:)=BPtemp;
    GOE(i,:)=GOEtemp;
    other(i,:)=othertemp;

end

% Format data:
for i=1:length(varvec)
    if min([CSS(i,2),BP(i,2),rep(i,2),other(i,2),GOE(i,2)])==min(CSS(i,2),BP(i,2)) && isnan(rep(i,2))
        rep(i,2)=0;
        rep(i,1)=i;
    end
end
if sum(~isnan(CSS),'all')~=0
    xvec=varvec(CSS(~isnan(CSS(:,1)),1));
    yvec=CSS(~isnan(CSS(:,2)),2);
else
    xvec=NaN(length(CSS(:,1)),1);
    yvec=NaN(length(CSS(:,1)),1);
end
if sum(~isnan(rep),'all')~=0
    xrepvec=varvec(rep(~isnan(rep(:,1)),1));
    yrepvec=rep(~isnan(rep(:,2)),2);
else
    xrepvec=NaN(length(rep(:,1)),1);
    yrepvec=NaN(length(rep(:,1)),1);
end
if sum(~isnan(BP),'all')~=0
    xBPvec=varvec(BP(~isnan(BP(:,1)),1));
    yBPvec=BP(~isnan(BP(:,2)),2);
else
    xBPvec=NaN(length(BP(:,1)),1);
    yBPvec=NaN(length(BP(:,1)),1);
end
if sum(~isnan(GOE),'all')~=0
    xGOEvec=varvec(GOE(~isnan(GOE(:,1)),1));
    yGOEvec=GOE(~isnan(GOE(:,2)),2);
else
    xGOEvec=NaN(length(GOE(:,1)),1);
    yGOEvec=NaN(length(GOE(:,1)),1);
end
if sum(~isnan(other),'all')~=0
    xothervec=varvec(other(~isnan(other(:,1)),1));
    yothervec=other(~isnan(other(:,2)),2);
else
    xothervec=NaN(length(other(:,1)),1);
    yothervec=NaN(length(other(:,1)),1);
end
xrepvec1=NaN(length(yrepvec),1);
xrepvec2=NaN(length(yrepvec),1);
yrepvec1=NaN(length(yrepvec),1);
yrepvec2=NaN(length(yrepvec),1);
for i=1:length(yrepvec)
    if yrepvec(i)>0.1
        yrepvec1(i)=yrepvec(i);
        xrepvec1(i)=xrepvec(i);
    else
        yrepvec2(i)=yrepvec(i);
        xrepvec2(i)=xrepvec(i);
    end
end
for i=1:length(varvec)
    if min([CSS(i,2),BP(i,2),rep(i,2),other(i,2),GOE(i,2)])==min(CSS(i,2),BP(i,2)) && ~isnan(rep(i,2))
        yrepvec2(i)=0;
        xrepvec2(i)=varvec(i);
    end
end
for i=1:length(varvec)
    if xrepvec1(i)==xrepvec2(i) && yrepvec1(i)==yrepvec2(i)
        xrepvec1(i)=NaN;
        yrepvec1(i)=NaN;
    end
end
xrepvec1(isnan(xrepvec1))=[];
xrepvec2(isnan(xrepvec2))=[];
yrepvec1(isnan(yrepvec1))=[];
yrepvec2(isnan(yrepvec2))=[];

% Create the plot:
subplot(3,3,3)
plot(xGOEvec,yGOEvec,'x','color',paleblue,'markersize',8,'linewidth',2)
axis square
hold on
plot(xBPvec,yBPvec,'--','color',purple,'markersize',8,'linewidth',3)
hold on
plot(xrepvec1,yrepvec1,':','color',orange,'markersize',8,'linewidth',3)
hold on
plot(xrepvec2,yrepvec2,':','color',orange,'markersize',8,'linewidth',3)
hold on
plot(xvec,yvec,'-','color',green,'markersize',8,'linewidth',3)
hold on
plot(xothervec,yothervec,'x','color',midblue,'markersize',8,'linewidth',2)
set(gca,'xtick',[0,1],'fontsize',14)
set(gca,'ytick',[0,1],'fontsize',14)
ylim([0,ymax])
text(0.05,0.9,'C','fontsize',16)
xlabel(xaxislabel,'interpreter','latex');
axis square

%% Whole-life trade-offs, effect of beta

% Define parameters:
a0=5;
c1a=0.1;
c2a=2;
b0=0.1;
c1b=0.1;
c2b=2;
f=1;
h=0.75;
alpha=1;
delta=0.25;
q=1;
gamma=0;
ymax=1;
xaxislabel='$\beta$';
version=5;
varvec=linspace(0.01,1,100);

% Set up vectors to use later:
CSS=NaN(length(varvec),2);
rep=NaN(length(varvec),2);
BP=NaN(length(varvec),2);
GOE=NaN(length(varvec),2);
other=NaN(length(varvec),2);
disp("working on plot 4")

% For each value of beta, we determine the singular strategies and their
% stability:
parfor i=1:length(varvec)
    
    % Set up vectors:
    CSStemp=NaN(1,2);
    reptemp=NaN(1,2);
    BPtemp=NaN(1,2);
    GOEtemp=NaN(1,2);
    othertemp=NaN(1,2);
    
    % For different combinations of trade-offs, we find the singular
    % strategies and their stability:
    beta=varvec(i);
    [CSSvec,repvec,BPvec,GOEvec,othervec]=find_singstrats_function(a0,c1a,c2a,b0,c1b,c2b,f,h,beta,alpha,delta,q,gamma,version,ymax);

    % Format output:
    if ~isempty(CSSvec)
        CSStemp=[i,CSSvec(1)];
    end
    if ~isempty(repvec)
        reptemp=[i,repvec(1)];
    end
    if ~isempty(BPvec)
        BPtemp=[i,BPvec(1)];
    end
    if ~isempty(GOEvec)
        GOEtemp=[i,GOEvec(1)];
    end
    if ~isempty(othervec)
        othertemp=[i,othervec(1)];
    end
    
    % It is very rare to see more than one singular strategy of the same
    % type. In these cases, the following messages appear:
    if length(CSSvec)>1
        disp("More than one CSS when i=" +i)
    end
    if length(repvec)>1
        disp("More than one repeller when i=" +i)
    end
    if length(BPvec)>1
        disp("More than one BP when i=" +i)
    end
    if length(GOEvec)>1
        disp("More than one GOE when i=" +i)
    end
    if length(othervec)>1
        disp("More than one 'other' when i=" +i)
    end
    
    % Create output:
    CSS(i,:)=CSStemp;
    rep(i,:)=reptemp;
    BP(i,:)=BPtemp;
    GOE(i,:)=GOEtemp;
    other(i,:)=othertemp;

end

% Format data:
for i=1:length(varvec)
    if min([CSS(i,2),BP(i,2),rep(i,2),other(i,2),GOE(i,2)])==min(CSS(i,2),BP(i,2)) && isnan(rep(i,2))
        rep(i,2)=0;
        rep(i,1)=i;
    end
end
if sum(~isnan(CSS),'all')~=0
    xvec=varvec(CSS(~isnan(CSS(:,1)),1));
    yvec=CSS(~isnan(CSS(:,2)),2);
else
    xvec=NaN(length(CSS(:,1)),1);
    yvec=NaN(length(CSS(:,1)),1);
end
if sum(~isnan(rep),'all')~=0
    xrepvec=varvec(rep(~isnan(rep(:,1)),1));
    yrepvec=rep(~isnan(rep(:,2)),2);
else
    xrepvec=NaN(length(rep(:,1)),1);
    yrepvec=NaN(length(rep(:,1)),1);
end
if sum(~isnan(BP),'all')~=0
    xBPvec=varvec(BP(~isnan(BP(:,1)),1));
    yBPvec=BP(~isnan(BP(:,2)),2);
else
    xBPvec=NaN(length(BP(:,1)),1);
    yBPvec=NaN(length(BP(:,1)),1);
end
if sum(~isnan(GOE),'all')~=0
    xGOEvec=varvec(GOE(~isnan(GOE(:,1)),1));
    yGOEvec=GOE(~isnan(GOE(:,2)),2);
else
    xGOEvec=NaN(length(GOE(:,1)),1);
    yGOEvec=NaN(length(GOE(:,1)),1);
end
if sum(~isnan(other),'all')~=0
    xothervec=varvec(other(~isnan(other(:,1)),1));
    yothervec=other(~isnan(other(:,2)),2);
else
    xothervec=NaN(length(other(:,1)),1);
    yothervec=NaN(length(other(:,1)),1);
end
xrepvec1=NaN(length(yrepvec),1);
xrepvec2=NaN(length(yrepvec),1);
yrepvec1=NaN(length(yrepvec),1);
yrepvec2=NaN(length(yrepvec),1);
for i=1:length(yrepvec)
    if yrepvec(i)>0.1
        yrepvec1(i)=yrepvec(i);
        xrepvec1(i)=xrepvec(i);
    else
        yrepvec2(i)=yrepvec(i);
        xrepvec2(i)=xrepvec(i);
    end
end
for i=1:length(varvec)
    if min([CSS(i,2),BP(i,2),rep(i,2),other(i,2),GOE(i,2)])==min(CSS(i,2),BP(i,2)) && ~isnan(rep(i,2))
        yrepvec2(i)=0;
        xrepvec2(i)=varvec(i);
    end
end
for i=1:length(varvec)
    if xrepvec1(i)==xrepvec2(i) && yrepvec1(i)==yrepvec2(i)
        xrepvec1(i)=NaN;
        yrepvec1(i)=NaN;
    end
end
xrepvec1(isnan(xrepvec1))=[];
xrepvec2(isnan(xrepvec2))=[];
yrepvec1(isnan(yrepvec1))=[];
yrepvec2(isnan(yrepvec2))=[];

% Create the plot:
subplot(3,3,4)
plot(xGOEvec,yGOEvec,'x','color',paleblue,'markersize',8,'linewidth',2)
axis square
hold on
plot(xBPvec,yBPvec,'--','color',purple,'markersize',8,'linewidth',3)
hold on
plot(xrepvec1,yrepvec1,':','color',orange,'markersize',8,'linewidth',3)
hold on
plot(xrepvec2,yrepvec2,':','color',orange,'markersize',8,'linewidth',3)
hold on
plot(xvec,yvec,'-','color',green,'markersize',8,'linewidth',3)
hold on
plot(xothervec,yothervec,'x','color',midblue,'markersize',8,'linewidth',2)
set(gca,'xtick',[0,1],'fontsize',14)
set(gca,'ytick',[0,1],'fontsize',14)
ylim([0,ymax])
text(0.05,0.9,'D','fontsize',16)
xlabel(xaxislabel,'interpreter','latex');
axis square
ylabel('Singular values of rate of onset of resistance, $\zeta$','fontsize',14,'interpreter','latex')

%% Whole-life trade-offs, effect of 1-f

% Define parameters:
a0=5;
c1a=0.1;
c2a=2;
b0=0.2;
c1b=0.1;
c2b=2;
beta=0.5;
h=0.75;
alpha=0;
delta=0.25;
q=1;
gamma=0;
ymax=1;
xaxislabel='$1-f$';
version=5;
varvec=linspace(0.01,1,100);

% Set up vectors to use later:
CSS=NaN(length(varvec),2);
rep=NaN(length(varvec),2);
BP=NaN(length(varvec),2);
GOE=NaN(length(varvec),2);
other=NaN(length(varvec),2);
disp("working on plot 5")

% For each value of 1-f, we determine the singular strategies and their
% stability:
parfor i=1:length(varvec)
    
    % Set up vectors:
    CSStemp=NaN(1,2);
    reptemp=NaN(1,2);
    BPtemp=NaN(1,2);
    GOEtemp=NaN(1,2);
    othertemp=NaN(1,2);
    
    % For different combinations of trade-offs, we find the singular
    % strategies and their stability:
    f=1-varvec(i);
    [CSSvec,repvec,BPvec,GOEvec,othervec]=find_singstrats_function(a0,c1a,c2a,b0,c1b,c2b,f,h,beta,alpha,delta,q,gamma,version,ymax);

    % Format output:
    if ~isempty(CSSvec)
        CSStemp=[i,CSSvec(1)];
    end
    if ~isempty(repvec)
        reptemp=[i,repvec(1)];
    end
    if ~isempty(BPvec)
        BPtemp=[i,BPvec(1)];
    end
    if ~isempty(GOEvec)
        GOEtemp=[i,GOEvec(1)];
    end
    if ~isempty(othervec)
        othertemp=[i,othervec(1)];
    end
    
    % It is very rare to see more than one singular strategy of the same
    % type. In these cases, the following messages appear:
    if length(CSSvec)>1
        disp("More than one CSS when i=" +i)
    end
    if length(repvec)>1
        disp("More than one repeller when i=" +i)
    end
    if length(BPvec)>1
        disp("More than one BP when i=" +i)
    end
    if length(GOEvec)>1
        disp("More than one GOE when i=" +i)
    end
    if length(othervec)>1
        disp("More than one 'other' when i=" +i)
    end
    
    % Create output:
    CSS(i,:)=CSStemp;
    rep(i,:)=reptemp;
    BP(i,:)=BPtemp;
    GOE(i,:)=GOEtemp;
    other(i,:)=othertemp;

end

% Format data:
for i=1:length(varvec)
    if min([CSS(i,2),BP(i,2),rep(i,2),other(i,2),GOE(i,2)])==min(CSS(i,2),BP(i,2)) && isnan(rep(i,2))
        rep(i,2)=0;
        rep(i,1)=i;
    end
end
if sum(~isnan(CSS),'all')~=0
    xvec=varvec(CSS(~isnan(CSS(:,1)),1));
    yvec=CSS(~isnan(CSS(:,2)),2);
else
    xvec=NaN(length(CSS(:,1)),1);
    yvec=NaN(length(CSS(:,1)),1);
end
if sum(~isnan(rep),'all')~=0
    xrepvec=varvec(rep(~isnan(rep(:,1)),1));
    yrepvec=rep(~isnan(rep(:,2)),2);
else
    xrepvec=NaN(length(rep(:,1)),1);
    yrepvec=NaN(length(rep(:,1)),1);
end
if sum(~isnan(BP),'all')~=0
    xBPvec=varvec(BP(~isnan(BP(:,1)),1));
    yBPvec=BP(~isnan(BP(:,2)),2);
else
    xBPvec=NaN(length(BP(:,1)),1);
    yBPvec=NaN(length(BP(:,1)),1);
end
if sum(~isnan(GOE),'all')~=0
    xGOEvec=varvec(GOE(~isnan(GOE(:,1)),1));
    yGOEvec=GOE(~isnan(GOE(:,2)),2);
else
    xGOEvec=NaN(length(GOE(:,1)),1);
    yGOEvec=NaN(length(GOE(:,1)),1);
end
if sum(~isnan(other),'all')~=0
    xothervec=varvec(other(~isnan(other(:,1)),1));
    yothervec=other(~isnan(other(:,2)),2);
else
    xothervec=NaN(length(other(:,1)),1);
    yothervec=NaN(length(other(:,1)),1);
end
xrepvec1=NaN(length(yrepvec),1);
xrepvec2=NaN(length(yrepvec),1);
yrepvec1=NaN(length(yrepvec),1);
yrepvec2=NaN(length(yrepvec),1);
for i=1:length(yrepvec)
    if yrepvec(i)>0.27
        yrepvec1(i)=yrepvec(i);
        xrepvec1(i)=xrepvec(i);
    else
        yrepvec2(i)=yrepvec(i);
        xrepvec2(i)=xrepvec(i);
    end
end
for i=1:length(varvec)
    if min([CSS(i,2),BP(i,2),rep(i,2),other(i,2),GOE(i,2)])==min(CSS(i,2),BP(i,2)) && ~isnan(rep(i,2))
        yrepvec2(i)=0;
        xrepvec2(i)=varvec(i);
    end
end
for i=1:length(varvec)
    if xrepvec1(i)==xrepvec2(i) && yrepvec1(i)==yrepvec2(i)
        xrepvec1(i)=NaN;
        yrepvec1(i)=NaN;
    end
end
xrepvec1(isnan(xrepvec1))=[];
xrepvec2(isnan(xrepvec2))=[];
yrepvec1(isnan(yrepvec1))=[];
yrepvec2(isnan(yrepvec2))=[];

% Create the plot:
subplot(3,3,5)
plot(xGOEvec,yGOEvec,'x','color',paleblue,'markersize',8,'linewidth',2)
axis square
hold on
plot(xBPvec,yBPvec,'--','color',purple,'markersize',8,'linewidth',3)
hold on
plot(xrepvec1,yrepvec1,':','color',orange,'markersize',8,'linewidth',3)
hold on
plot(xrepvec2,yrepvec2,':','color',orange,'markersize',8,'linewidth',3)
hold on
plot(xvec,yvec,'-','color',green,'markersize',8,'linewidth',3)
hold on
plot(xothervec,yothervec,'x','color',midblue,'markersize',8,'linewidth',2)
set(gca,'xtick',[0,1],'fontsize',14)
set(gca,'ytick',[0,1],'fontsize',14)
ylim([0,ymax])
text(0.05,0.9,'E','fontsize',16)
xlabel(xaxislabel,'interpreter','latex');
axis square

%% Whole-life trade-offs, effect of alpha

% Define parameters:
a0=5;
c1a=0.2;
c2a=10;
b0=0.1;
c1b=0.2;
c2b=10;
h=0.75;
f=1;
delta=0.25;
beta=0.5;
q=1;
gamma=0;
ymax=1;
xaxislabel='$\alpha$';
version=5;
varvec=linspace(0.04,4,100);

% Set up vectors to use later:
CSS=NaN(length(varvec),2);
rep=NaN(length(varvec),2);
BP=NaN(length(varvec),2);
GOE=NaN(length(varvec),2);
other=NaN(length(varvec),2);
disp("working on plot 6")

% For each value of alpha, we determine the singular strategies and their
% stability:
parfor i=1:length(varvec)
    
    % Set up vectors:
    CSStemp=NaN(1,2);
    reptemp=NaN(1,2);
    BPtemp=NaN(1,2);
    GOEtemp=NaN(1,2);
    othertemp=NaN(1,2);
    
    % For different combinations of trade-offs, we find the singular
    % strategies and their stability:
    alpha=varvec(i);
    [CSSvec,repvec,BPvec,GOEvec,othervec]=find_singstrats_function(a0,c1a,c2a,b0,c1b,c2b,f,h,beta,alpha,delta,q,gamma,version,ymax);

    % Format output:
    if ~isempty(CSSvec)
        CSStemp=[i,CSSvec(1)];
    end
    if ~isempty(repvec)
        reptemp=[i,repvec(1)];
    end
    if ~isempty(BPvec)
        BPtemp=[i,BPvec(1)];
    end
    if ~isempty(GOEvec)
        GOEtemp=[i,GOEvec(1)];
    end
    if ~isempty(othervec)
        othertemp=[i,othervec(1)];
    end
    
    % It is very rare to see more than one singular strategy of the same
    % type. In these cases, the following messages appear:
    if length(CSSvec)>1
        disp("More than one CSS when i=" +i)
    end
    if length(repvec)>1
        disp("More than one repeller when i=" +i)
    end
    if length(BPvec)>1
        disp("More than one BP when i=" +i)
    end
    if length(GOEvec)>1
        disp("More than one GOE when i=" +i)
    end
    if length(othervec)>1
        disp("More than one 'other' when i=" +i)
    end
    
    % Create output:
    CSS(i,:)=CSStemp;
    rep(i,:)=reptemp;
    BP(i,:)=BPtemp;
    GOE(i,:)=GOEtemp;
    other(i,:)=othertemp;

end

% Format data:
for i=1:length(varvec)
    if min([CSS(i,2),BP(i,2),rep(i,2),other(i,2),GOE(i,2)])==min(CSS(i,2),BP(i,2)) && isnan(rep(i,2))
        rep(i,2)=0;
        rep(i,1)=i;
    end
end
if sum(~isnan(CSS),'all')~=0
    xvec=varvec(CSS(~isnan(CSS(:,1)),1));
    yvec=CSS(~isnan(CSS(:,2)),2);
else
    xvec=NaN(length(CSS(:,1)),1);
    yvec=NaN(length(CSS(:,1)),1);
end
if sum(~isnan(rep),'all')~=0
    xrepvec=varvec(rep(~isnan(rep(:,1)),1));
    yrepvec=rep(~isnan(rep(:,2)),2);
else
    xrepvec=NaN(length(rep(:,1)),1);
    yrepvec=NaN(length(rep(:,1)),1);
end
if sum(~isnan(BP),'all')~=0
    xBPvec=varvec(BP(~isnan(BP(:,1)),1));
    yBPvec=BP(~isnan(BP(:,2)),2);
else
    xBPvec=NaN(length(BP(:,1)),1);
    yBPvec=NaN(length(BP(:,1)),1);
end
if sum(~isnan(GOE),'all')~=0
    xGOEvec=varvec(GOE(~isnan(GOE(:,1)),1));
    yGOEvec=GOE(~isnan(GOE(:,2)),2);
else
    xGOEvec=NaN(length(GOE(:,1)),1);
    yGOEvec=NaN(length(GOE(:,1)),1);
end
if sum(~isnan(other),'all')~=0
    xothervec=varvec(other(~isnan(other(:,1)),1));
    yothervec=other(~isnan(other(:,2)),2);
else
    xothervec=NaN(length(other(:,1)),1);
    yothervec=NaN(length(other(:,1)),1);
end
xrepvec1=NaN(length(yrepvec),1);
xrepvec2=NaN(length(yrepvec),1);
yrepvec1=NaN(length(yrepvec),1);
yrepvec2=NaN(length(yrepvec),1);
for i=1:length(yrepvec)
    if yrepvec(i)>0.2-xrepvec(i)/15
        yrepvec1(i)=yrepvec(i);
        xrepvec1(i)=xrepvec(i);
    else
        yrepvec2(i)=yrepvec(i);
        xrepvec2(i)=xrepvec(i);
    end
end
for i=1:length(varvec)
    if min([CSS(i,2),BP(i,2),rep(i,2),other(i,2),GOE(i,2)])==min(CSS(i,2),BP(i,2)) && ~isnan(rep(i,2))
        yrepvec2(i)=0;
        xrepvec2(i)=varvec(i);
    end
end
for i=1:length(varvec)
    if xrepvec1(i)==xrepvec2(i) && yrepvec1(i)==yrepvec2(i)
        xrepvec1(i)=NaN;
        yrepvec1(i)=NaN;
    end
end
xrepvec1(isnan(xrepvec1))=[];
xrepvec2(isnan(xrepvec2))=[];
yrepvec1(isnan(yrepvec1))=[];
yrepvec2(isnan(yrepvec2))=[];

% Create the plot:
subplot(3,3,6)
plot(xGOEvec,yGOEvec,'x','color',paleblue,'markersize',8,'linewidth',2)
axis square
hold on
plot(xBPvec,yBPvec,'--','color',purple,'markersize',8,'linewidth',3)
hold on
plot(xrepvec1,yrepvec1,':','color',orange,'markersize',8,'linewidth',3)
hold on
plot(xrepvec2,yrepvec2,':','color',orange,'markersize',8,'linewidth',3)
hold on
plot(xvec,yvec,'-','color',green,'markersize',8,'linewidth',3)
hold on
plot(xothervec,yothervec,'x','color',midblue,'markersize',8,'linewidth',2)
set(gca,'xtick',[0,4],'fontsize',14)
set(gca,'ytick',[0,1],'fontsize',14)
ylim([0,ymax])
text(0.2,0.9,'F','fontsize',16)
xlabel(xaxislabel,'interpreter','latex');
axis square

%% Whole-life trade-offs, effect of b0

% Define parameters:
a0=5;
c1a=0.1;
c2a=2;
c1b=0.1;
c2b=2;
h=0.75;
alpha=0;
f=0.5;
delta=0.25;
beta=0.5;
q=1;
gamma=0;
ymax=1;
xaxislabel='$b_0$';
version=5;
varvec=linspace(0.005,0.5,100);

% Set up vectors to use later:
CSS=NaN(length(varvec),2);
rep=NaN(length(varvec),2);
BP=NaN(length(varvec),2);
GOE=NaN(length(varvec),2);
other=NaN(length(varvec),2);
disp("working on plot 7")

% For each value of b0, we determine the singular strategies and their
% stability:
parfor i=1:length(varvec)
    
    % Set up vectors:
    CSStemp=NaN(1,2);
    reptemp=NaN(1,2);
    BPtemp=NaN(1,2);
    GOEtemp=NaN(1,2);
    othertemp=NaN(1,2);
    
    % For different combinations of trade-offs, we find the singular
    % strategies and their stability:
    b0=varvec(i);
    [CSSvec,repvec,BPvec,GOEvec,othervec]=find_singstrats_function(a0,c1a,c2a,b0,c1b,c2b,f,h,beta,alpha,delta,q,gamma,version,ymax);

    % Format output:
    if ~isempty(CSSvec)
        CSStemp=[i,CSSvec(1)];
    end
    if ~isempty(repvec)
        reptemp=[i,repvec(1)];
    end
    if ~isempty(BPvec)
        BPtemp=[i,BPvec(1)];
    end
    if ~isempty(GOEvec)
        GOEtemp=[i,GOEvec(1)];
    end
    if ~isempty(othervec)
        othertemp=[i,othervec(1)];
    end
    
    % It is very rare to see more than one singular strategy of the same
    % type. In these cases, the following messages appear:
    if length(CSSvec)>1
        disp("More than one CSS when i=" +i)
    end
    if length(repvec)>1
        disp("More than one repeller when i=" +i)
    end
    if length(BPvec)>1
        disp("More than one BP when i=" +i)
    end
    if length(GOEvec)>1
        disp("More than one GOE when i=" +i)
    end
    if length(othervec)>1
        disp("More than one 'other' when i=" +i)
    end
    
    % Create output:
    CSS(i,:)=CSStemp;
    rep(i,:)=reptemp;
    BP(i,:)=BPtemp;
    GOE(i,:)=GOEtemp;
    other(i,:)=othertemp;

end

% Format data:
for i=1:length(varvec)
    if min([CSS(i,2),BP(i,2),rep(i,2),other(i,2),GOE(i,2)])==min(CSS(i,2),BP(i,2)) && isnan(rep(i,2))
        rep(i,2)=0;
        rep(i,1)=i;
    end
end
if sum(~isnan(CSS),'all')~=0
    xvec=varvec(CSS(~isnan(CSS(:,1)),1));
    yvec=CSS(~isnan(CSS(:,2)),2);
else
    xvec=NaN(length(CSS(:,1)),1);
    yvec=NaN(length(CSS(:,1)),1);
end
if sum(~isnan(rep),'all')~=0
    xrepvec=varvec(rep(~isnan(rep(:,1)),1));
    yrepvec=rep(~isnan(rep(:,2)),2);
else
    xrepvec=NaN(length(rep(:,1)),1);
    yrepvec=NaN(length(rep(:,1)),1);
end
if sum(~isnan(BP),'all')~=0
    xBPvec=varvec(BP(~isnan(BP(:,1)),1));
    yBPvec=BP(~isnan(BP(:,2)),2);
else
    xBPvec=NaN(length(BP(:,1)),1);
    yBPvec=NaN(length(BP(:,1)),1);
end
if sum(~isnan(GOE),'all')~=0
    xGOEvec=varvec(GOE(~isnan(GOE(:,1)),1));
    yGOEvec=GOE(~isnan(GOE(:,2)),2);
else
    xGOEvec=NaN(length(GOE(:,1)),1);
    yGOEvec=NaN(length(GOE(:,1)),1);
end
if sum(~isnan(other),'all')~=0
    xothervec=varvec(other(~isnan(other(:,1)),1));
    yothervec=other(~isnan(other(:,2)),2);
else
    xothervec=NaN(length(other(:,1)),1);
    yothervec=NaN(length(other(:,1)),1);
end
xrepvec1=NaN(length(yrepvec),1);
xrepvec2=NaN(length(yrepvec),1);
yrepvec1=NaN(length(yrepvec),1);
yrepvec2=NaN(length(yrepvec),1);
for i=1:length(yrepvec)
    if yrepvec(i)>0.1
        yrepvec1(i)=yrepvec(i);
        xrepvec1(i)=xrepvec(i);
    else
        yrepvec2(i)=yrepvec(i);
        xrepvec2(i)=xrepvec(i);
    end
end
for i=1:length(varvec)
    if min([CSS(i,2),BP(i,2),rep(i,2),other(i,2),GOE(i,2)])==min(CSS(i,2),BP(i,2)) && ~isnan(rep(i,2))
        yrepvec2(i)=0;
        xrepvec2(i)=varvec(i);
    end
end
for i=1:length(varvec)
    if xrepvec1(i)==xrepvec2(i) && yrepvec1(i)==yrepvec2(i)
        xrepvec1(i)=NaN;
        yrepvec1(i)=NaN;
    end
end
xrepvec1(isnan(xrepvec1))=[];
xrepvec2(isnan(xrepvec2))=[];
yrepvec1(isnan(yrepvec1))=[];
yrepvec2(isnan(yrepvec2))=[];

% Create the plot:
subplot(3,3,7)
plot(xGOEvec,yGOEvec,'x','color',paleblue,'markersize',8,'linewidth',2)
axis square
hold on
plot(xBPvec,yBPvec,'--','color',purple,'markersize',8,'linewidth',3)
hold on
plot(xrepvec1,yrepvec1,':','color',orange,'markersize',8,'linewidth',3)
hold on
plot(xrepvec2,yrepvec2,':','color',orange,'markersize',8,'linewidth',3)
hold on
plot(xvec,yvec,'-','color',green,'markersize',8,'linewidth',3)
hold on
plot(xothervec,yothervec,'x','color',midblue,'markersize',8,'linewidth',2)
set(gca,'xtick',[0,0.5],'fontsize',14)
set(gca,'ytick',[0,1],'fontsize',14)
ylim([0,ymax])
text(0.025,0.9,'G','fontsize',16)
xlabel(xaxislabel,'interpreter','latex');
axis square

%% Whole-life trade-offs, effect of c1a

% Define parameters:
a0=5;
c2a=2;
b0=0.1;
c1b=0.1;
c2b=2;
h=0.75;
alpha=1;
f=1;
delta=0.25;
beta=0.5;
q=1;
gamma=0;
ymax=1;
xaxislabel='$c_1^a$';
version=5;
varvec=linspace(0.01,1,100);

% Set up vectors to use later:
CSS=NaN(length(varvec),2);
rep=NaN(length(varvec),2);
BP=NaN(length(varvec),2);
GOE=NaN(length(varvec),2);
other=NaN(length(varvec),2);
disp("working on plot 8")

% For each value of c1a, we determine the singular strategies and their
% stability:
parfor i=1:length(varvec)
    
    % Set up vectors:
    CSStemp=NaN(1,2);
    reptemp=NaN(1,2);
    BPtemp=NaN(1,2);
    GOEtemp=NaN(1,2);
    othertemp=NaN(1,2);
    
    % For different combinations of trade-offs, we find the singular
    % strategies and their stability:
    c1a=varvec(i);
    [CSSvec,repvec,BPvec,GOEvec,othervec]=find_singstrats_function(a0,c1a,c2a,b0,c1b,c2b,f,h,beta,alpha,delta,q,gamma,version,ymax);

    % Format output:
    if ~isempty(CSSvec)
        CSStemp=[i,CSSvec(1)];
    end
    if ~isempty(repvec)
        reptemp=[i,repvec(1)];
    end
    if ~isempty(BPvec)
        BPtemp=[i,BPvec(1)];
    end
    if ~isempty(GOEvec)
        GOEtemp=[i,GOEvec(1)];
    end
    if ~isempty(othervec)
        othertemp=[i,othervec(1)];
    end
    
    % It is very rare to see more than one singular strategy of the same
    % type. In these cases, the following messages appear:
    if length(CSSvec)>1
        disp("More than one CSS when i=" +i)
    end
    if length(repvec)>1
        disp("More than one repeller when i=" +i)
    end
    if length(BPvec)>1
        disp("More than one BP when i=" +i)
    end
    if length(GOEvec)>1
        disp("More than one GOE when i=" +i)
    end
    if length(othervec)>1
        disp("More than one 'other' when i=" +i)
    end
    
    % Create output:
    CSS(i,:)=CSStemp;
    rep(i,:)=reptemp;
    BP(i,:)=BPtemp;
    GOE(i,:)=GOEtemp;
    other(i,:)=othertemp;

end

% Format data:
for i=1:length(varvec)
    if min([CSS(i,2),BP(i,2),rep(i,2),other(i,2),GOE(i,2)])==min(CSS(i,2),BP(i,2)) && isnan(rep(i,2))
        rep(i,2)=0;
        rep(i,1)=i;
    end
end
if sum(~isnan(CSS),'all')~=0
    xvec=varvec(CSS(~isnan(CSS(:,1)),1));
    yvec=CSS(~isnan(CSS(:,2)),2);
else
    xvec=NaN(length(CSS(:,1)),1);
    yvec=NaN(length(CSS(:,1)),1);
end
if sum(~isnan(rep),'all')~=0
    xrepvec=varvec(rep(~isnan(rep(:,1)),1));
    yrepvec=rep(~isnan(rep(:,2)),2);
else
    xrepvec=NaN(length(rep(:,1)),1);
    yrepvec=NaN(length(rep(:,1)),1);
end
if sum(~isnan(BP),'all')~=0
    xBPvec=varvec(BP(~isnan(BP(:,1)),1));
    yBPvec=BP(~isnan(BP(:,2)),2);
else
    xBPvec=NaN(length(BP(:,1)),1);
    yBPvec=NaN(length(BP(:,1)),1);
end
if sum(~isnan(GOE),'all')~=0
    xGOEvec=varvec(GOE(~isnan(GOE(:,1)),1));
    yGOEvec=GOE(~isnan(GOE(:,2)),2);
else
    xGOEvec=NaN(length(GOE(:,1)),1);
    yGOEvec=NaN(length(GOE(:,1)),1);
end
if sum(~isnan(other),'all')~=0
    xothervec=varvec(other(~isnan(other(:,1)),1));
    yothervec=other(~isnan(other(:,2)),2);
else
    xothervec=NaN(length(other(:,1)),1);
    yothervec=NaN(length(other(:,1)),1);
end
xrepvec1=NaN(length(yrepvec),1);
xrepvec2=NaN(length(yrepvec),1);
yrepvec1=NaN(length(yrepvec),1);
yrepvec2=NaN(length(yrepvec),1);
for i=1:length(yrepvec)
    if yrepvec(i)>0.1
        yrepvec1(i)=yrepvec(i);
        xrepvec1(i)=xrepvec(i);
    else
        yrepvec2(i)=yrepvec(i);
        xrepvec2(i)=xrepvec(i);
    end
end
for i=1:length(varvec)
    if min([CSS(i,2),BP(i,2),rep(i,2),other(i,2),GOE(i,2)])==min(CSS(i,2),BP(i,2)) && ~isnan(rep(i,2))
        yrepvec2(i)=0;
        xrepvec2(i)=varvec(i);
    end
end
for i=1:length(varvec)
    if xrepvec1(i)==xrepvec2(i) && yrepvec1(i)==yrepvec2(i)
        xrepvec1(i)=NaN;
        yrepvec1(i)=NaN;
    end
end
xrepvec1(isnan(xrepvec1))=[];
xrepvec2(isnan(xrepvec2))=[];
yrepvec1(isnan(yrepvec1))=[];
yrepvec2(isnan(yrepvec2))=[];

% Create the plot:
subplot(3,3,8)
plot(xGOEvec,yGOEvec,'x','color',paleblue,'markersize',8,'linewidth',2)
axis square
hold on
plot(xBPvec,yBPvec,'--','color',purple,'markersize',8,'linewidth',3)
hold on
plot(xrepvec1,yrepvec1,':','color',orange,'markersize',8,'linewidth',3)
hold on
plot(xrepvec2,yrepvec2,':','color',orange,'markersize',8,'linewidth',3)
hold on
plot(xvec,yvec,'-','color',green,'markersize',8,'linewidth',3)
hold on
plot(xothervec,yothervec,'x','color',midblue,'markersize',8,'linewidth',2)
set(gca,'xtick',[0,1],'fontsize',14)
set(gca,'ytick',[0,1],'fontsize',14)
ylim([0,ymax])
text(0.025,0.9,'H','fontsize',16)
xlabel(xaxislabel,'interpreter','latex');
axis square

%% Whole-life trade-offs, effect of c1b

% Define parameters:
a0=5;
c2a=8;
b0=0.1;
c1a=0.1;
c2b=8;
h=0.75;
alpha=1;
f=1;
delta=0.25;
beta=0.5;
q=1;
gamma=0;
ymax=1;
xaxislabel='$c_1^b$';
version=6;
varvec=linspace(0.01,1,100);

% Set up vectors to use later:
CSS=NaN(length(varvec),2);
rep=NaN(length(varvec),2);
BP=NaN(length(varvec),2);
GOE=NaN(length(varvec),2);
other=NaN(length(varvec),2);
disp("working on plot 9")

% For each value of c1b, we determine the singular strategies and their
% stability:
parfor i=1:length(varvec)
    
    % Set up vectors:
    CSStemp=NaN(1,2);
    reptemp=NaN(1,2);
    BPtemp=NaN(1,2);
    GOEtemp=NaN(1,2);
    othertemp=NaN(1,2);
    
    % For different combinations of trade-offs, we find the singular
    % strategies and their stability:
    c1b=varvec(i);
    [CSSvec,repvec,BPvec,GOEvec,othervec]=find_singstrats_function(a0,c1a,c2a,b0,c1b,c2b,f,h,beta,alpha,delta,q,gamma,version,ymax);

    % Format output:
    if ~isempty(CSSvec)
        CSStemp=[i,CSSvec(1)];
    end
    if ~isempty(repvec)
        reptemp=[i,repvec(1)];
    end
    if ~isempty(BPvec)
        BPtemp=[i,BPvec(1)];
    end
    if ~isempty(GOEvec)
        GOEtemp=[i,GOEvec(1)];
    end
    if ~isempty(othervec)
        othertemp=[i,othervec(1)];
    end
    
    % It is very rare to see more than one singular strategy of the same
    % type. In these cases, the following messages appear:
    if length(CSSvec)>1
        disp("More than one CSS when i=" +i)
    end
    if length(repvec)>1
        disp("More than one repeller when i=" +i)
    end
    if length(BPvec)>1
        disp("More than one BP when i=" +i)
    end
    if length(GOEvec)>1
        disp("More than one GOE when i=" +i)
    end
    if length(othervec)>1
        disp("More than one 'other' when i=" +i)
    end
    
    % Create output:
    CSS(i,:)=CSStemp;
    rep(i,:)=reptemp;
    BP(i,:)=BPtemp;
    GOE(i,:)=GOEtemp;
    other(i,:)=othertemp;

end

% Format data:
for i=1:length(varvec)
    if min([CSS(i,2),BP(i,2),rep(i,2),other(i,2),GOE(i,2)])==min(CSS(i,2),BP(i,2)) && isnan(rep(i,2))
        rep(i,2)=0;
        rep(i,1)=i;
    end
end
if sum(~isnan(CSS),'all')~=0
    xvec=varvec(CSS(~isnan(CSS(:,1)),1));
    yvec=CSS(~isnan(CSS(:,2)),2);
else
    xvec=NaN(length(CSS(:,1)),1);
    yvec=NaN(length(CSS(:,1)),1);
end
if sum(~isnan(rep),'all')~=0
    xrepvec=varvec(rep(~isnan(rep(:,1)),1));
    yrepvec=rep(~isnan(rep(:,2)),2);
else
    xrepvec=NaN(length(rep(:,1)),1);
    yrepvec=NaN(length(rep(:,1)),1);
end
if sum(~isnan(BP),'all')~=0
    xBPvec=varvec(BP(~isnan(BP(:,1)),1));
    yBPvec=BP(~isnan(BP(:,2)),2);
else
    xBPvec=NaN(length(BP(:,1)),1);
    yBPvec=NaN(length(BP(:,1)),1);
end
if sum(~isnan(GOE),'all')~=0
    xGOEvec=varvec(GOE(~isnan(GOE(:,1)),1));
    yGOEvec=GOE(~isnan(GOE(:,2)),2);
else
    xGOEvec=NaN(length(GOE(:,1)),1);
    yGOEvec=NaN(length(GOE(:,1)),1);
end
if sum(~isnan(other),'all')~=0
    xothervec=varvec(other(~isnan(other(:,1)),1));
    yothervec=other(~isnan(other(:,2)),2);
else
    xothervec=NaN(length(other(:,1)),1);
    yothervec=NaN(length(other(:,1)),1);
end
xrepvec1=NaN(length(yrepvec),1);
xrepvec2=NaN(length(yrepvec),1);
yrepvec1=NaN(length(yrepvec),1);
yrepvec2=NaN(length(yrepvec),1);
for i=1:length(yrepvec)
    if yrepvec(i)>0.1
        yrepvec1(i)=yrepvec(i);
        xrepvec1(i)=xrepvec(i);
    else
        yrepvec2(i)=yrepvec(i);
        xrepvec2(i)=xrepvec(i);
    end
end
for i=1:length(varvec)
    if min([CSS(i,2),BP(i,2),rep(i,2),other(i,2),GOE(i,2)])==min(CSS(i,2),BP(i,2)) && ~isnan(rep(i,2))
        yrepvec2(i)=0;
        xrepvec2(i)=varvec(i);
    end
end
for i=1:length(varvec)
    if xrepvec1(i)==xrepvec2(i) && yrepvec1(i)==yrepvec2(i)
        xrepvec1(i)=NaN;
        yrepvec1(i)=NaN;
    end
end
xrepvec1(isnan(xrepvec1))=[];
xrepvec2(isnan(xrepvec2))=[];
yrepvec1(isnan(yrepvec1))=[];
yrepvec2(isnan(yrepvec2))=[];

% Create the plot:
subplot(3,3,9)
plot(xGOEvec,yGOEvec,'x','color',paleblue,'markersize',8,'linewidth',2)
axis square
hold on
plot(xBPvec,yBPvec,'--','color',purple,'markersize',8,'linewidth',3)
hold on
plot(xrepvec1,yrepvec1,':','color',orange,'markersize',8,'linewidth',3)
hold on
plot(xrepvec2,yrepvec2,':','color',orange,'markersize',8,'linewidth',3)
hold on
plot(xvec,yvec,'-','color',green,'markersize',8,'linewidth',3)
hold on
plot(xothervec,yothervec,'x','color',midblue,'markersize',8,'linewidth',2)
set(gca,'xtick',[0,1],'fontsize',14)
set(gca,'ytick',[0,1],'fontsize',14)
ylim([0,ymax])
text(0.025,0.9,'I','fontsize',16)
xlabel(xaxislabel,'interpreter','latex');
axis square
