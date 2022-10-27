% This code plots curves showing the effect of trade-off shape (c2a) on 
% singular strategies.

% Define parameters:
a0=5;
c1a=0.1;
c2a=2;
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
xaxislabel='c2a';

% We will vary the trade-off shape, c2a:
varvec=linspace(0.1,10,100);

% Set up vectors to use later:
CSS2=NaN(length(varvec),2);
rep2=NaN(length(varvec),2);
BP2=NaN(length(varvec),2);
GOE2=NaN(length(varvec),2);
other2=NaN(length(varvec),2);
CSS3=NaN(length(varvec),2);
rep3=NaN(length(varvec),2);
BP3=NaN(length(varvec),2);
GOE3=NaN(length(varvec),2);
other3=NaN(length(varvec),2);

% For each value of c2a, we determine the singular strategies and their
% stability:
parfor i=1:length(varvec)
    disp(i)
    
    % Set up vectors:
    CSS2temp=NaN(1,2);
    rep2temp=NaN(1,2);
    BP2temp=NaN(1,2);
    GOE2temp=NaN(1,2);
    other2temp=NaN(1,2);
    CSS3temp=NaN(1,2);
    rep3temp=NaN(1,2);
    BP3temp=NaN(1,2);
    GOE3temp=NaN(1,2);
    other3temp=NaN(1,2);
    
    % For different combinations of trade-offs, we find the singular
    % strategies and their stability:
    c2a=varvec(i);
    version=3;
    [CSSvec2,repvec2,BPvec2,GOEvec2,othervec2]=find_singstrats_function(a0,c1a,c2a,b0,c1b,c2b,f,h,beta,alpha,delta,q,gamma,version,ymax);
    version=5;
    [CSSvec3,repvec3,BPvec3,GOEvec3,othervec3]=find_singstrats_function(a0,c1a,c2a,b0,c1b,c2b,f,h,beta,alpha,delta,q,gamma,version,ymax);

    % Format output:
    if ~isempty(CSSvec2)
        CSS2temp=[i,CSSvec2(1)];
    end
    if ~isempty(repvec2)
        rep2temp=[i,repvec2(1)];
    end
    if ~isempty(BPvec2)
        BP2temp=[i,BPvec2(1)];
    end
    if ~isempty(GOEvec2)
        GOE2temp=[i,GOEvec2(1)];
    end
    if ~isempty(othervec2)
        other2temp=[i,othervec2(1)];
    end
    if ~isempty(CSSvec3)
        CSS3temp=[i,CSSvec3(1)];
    end
    if ~isempty(repvec3)
        rep3temp=[i,repvec3(1)];
    end
    if ~isempty(BPvec3)
        BP3temp=[i,BPvec3(1)];
    end
    if ~isempty(GOEvec3)
        GOE3temp=[i,GOEvec3(1)];
    end
    if ~isempty(othervec3)
        other3temp=[i,othervec3(1)];
    end
    
    % It is very rare to see more than one singular strategy of the same
    % type. In these cases, the following messages appear:
    if length(CSSvec2)>1
        disp("More than one CSS when version=2 and i=" +i)
    end
    if length(repvec2)>1
        disp("More than one repeller when version=2 and i=" +i)
    end
    if length(BPvec2)>1
        disp("More than one BP when version=2 and i=" +i)
    end
    if length(GOEvec2)>1
        disp("More than one GOE when version=2 and i=" +i)
    end
    if length(othervec2)>1
        disp("More than one 'other' when version=2 and i=" +i)
    end
    if length(CSSvec3)>1
        disp("More than one CSS when version=3 and i=" +i)
    end
    if length(repvec3)>1
        disp("More than one repeller when version=3 and i=" +i)
    end
    if length(BPvec3)>1
        disp("More than one BP when version=3 and i=" +i)
    end
    if length(GOEvec3)>1
        disp("More than one GOE when version=3 and i=" +i)
    end
    if length(othervec3)>1
        disp("More than one 'other' when version=3 and i=" +i)
    end
    
    % Create output:
    CSS2(i,:)=CSS2temp;
    rep2(i,:)=rep2temp;
    BP2(i,:)=BP2temp;
    GOE2(i,:)=GOE2temp;
    other2(i,:)=other2temp;
    CSS3(i,:)=CSS3temp;
    rep3(i,:)=rep3temp;
    BP3(i,:)=BP3temp;
    GOE3(i,:)=GOE3temp;
    other3(i,:)=other3temp;

end

%% Create the plot

% Colours for plotting:
green=1/255*[27,158,119];
orange=1/255*[217,95,2];
purple=1/255*[117,112,179];
paleblue=[1,0.6,1];
midblue=[0,0.6,0.1];

% Plot 1: Format data and create the plot
subplot(1,2,1)
if sum(~isnan(CSS2),'all')~=0
    xvec=varvec(CSS2(~isnan(CSS2(:,1)),1));
    yvec=CSS2(~isnan(CSS2(:,2)),2);
else
    xvec=NaN(length(CSS2(:,1)),1);
    yvec=NaN(length(CSS2(:,1)),1);
end
if sum(~isnan(rep2),'all')~=0
    xrepvec=varvec(rep2(~isnan(rep2(:,1)),1));
    yrepvec=rep2(~isnan(rep2(:,2)),2);
else
    xrepvec=NaN(length(rep2(:,1)),1);
    yrepvec=NaN(length(rep2(:,1)),1);
end
if sum(~isnan(BP2),'all')~=0
    xBPvec=varvec(BP2(~isnan(BP2(:,1)),1));
    yBPvec=BP2(~isnan(BP2(:,2)),2);
else
    xBPvec=NaN(length(BP2(:,1)),1);
    yBPvec=NaN(length(BP2(:,1)),1);
end
if sum(~isnan(GOE2),'all')~=0
    xGOEvec=varvec(GOE2(~isnan(GOE2(:,1)),1));
    yGOEvec=GOE2(~isnan(GOE2(:,2)),2);
else
    xGOEvec=NaN(length(GOE2(:,1)),1);
    yGOEvec=NaN(length(GOE2(:,1)),1);
end
if sum(~isnan(other2),'all')~=0
    xothervec=varvec(other2(~isnan(other2(:,1)),1));
    yothervec=other2(~isnan(other2(:,2)),2);
else
    xothervec=NaN(length(other2(:,1)),1);
    yothervec=NaN(length(other2(:,1)),1);
end
plot(xGOEvec,yGOEvec,'x','color',paleblue,'markersize',8,'linewidth',2)
axis square
hold on
plot(xBPvec,yBPvec,'--','color',purple,'markersize',8,'linewidth',3)
hold on
plot(xrepvec,yrepvec,':','color',orange,'markersize',8,'linewidth',3)
hold on
plot(xvec,yvec,'-','color',green,'markersize',8,'linewidth',3)
hold on
plot(xothervec,yothervec,'x','color',midblue,'markersize',8,'linewidth',2)
ylim([0,ymax])
set(gca,'xtick',[0,2,4,6,8,10],'fontsize',14)
set(gca,'ytick',[0,0.2,0.4,0.6,0.8,1],'fontsize',14)
ylabel('Singular strategies, $\zeta$', 'interpreter','latex','fontsize',20)
xlabel('Shape of trade-offs, $c_2^a$','interpreter','latex','fontsize',20)
title('Whole-life reproduction trade-off','fontsize',18,'interpreter','latex')
text(0.5,0.95,'A','fontsize',20)

% Plot 2: Format data and create the plot
subplot(1,2,2)
if sum(~isnan(CSS3),'all')~=0
    xvec=varvec(CSS3(~isnan(CSS3(:,1)),1));
    yvec=CSS3(~isnan(CSS3(:,2)),2);
else
    xvec=NaN(length(CSS3(:,1)),1);
    yvec=NaN(length(CSS3(:,1)),1);
end
if sum(~isnan(rep3),'all')~=0
    xrepvec=varvec(rep3(~isnan(rep3(:,1)),1));
    yrepvec=rep3(~isnan(rep3(:,2)),2);
else
    xrepvec=NaN(length(rep3(:,1)),1);
    yrepvec=NaN(length(rep3(:,1)),1);
end
if sum(~isnan(BP3),'all')~=0
    xBPvec=varvec(BP3(~isnan(BP3(:,1)),1));
    yBPvec=BP3(~isnan(BP3(:,2)),2);
else
    xBPvec=NaN(length(BP3(:,1)),1);
    yBPvec=NaN(length(BP3(:,1)),1);
end
if sum(~isnan(GOE3),'all')~=0
    xGOEvec=varvec(GOE3(~isnan(GOE3(:,1)),1));
    yGOEvec=GOE3(~isnan(GOE3(:,2)),2);
else
    xGOEvec=NaN(length(GOE3(:,1)),1);
    yGOEvec=NaN(length(GOE3(:,1)),1);
end
if sum(~isnan(other3),'all')~=0
    xothervec=varvec(other3(~isnan(other3(:,1)),1));
    yothervec=other3(~isnan(other3(:,2)),2);
else
    xothervec=NaN(length(other3(:,1)),1);
    yothervec=NaN(length(other3(:,1)),1);
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
set(gca,'xtick',[0,2,4,6,8,10],'fontsize',14)
set(gca,'ytick',[0,0.2,0.4,0.6,0.8,1],'fontsize',14)
ylim([0,ymax])
title('Pre-resistance reproduction trade-off','fontsize',18,'interpreter','latex')
text(0.5,0.95,'B','fontsize',20)
