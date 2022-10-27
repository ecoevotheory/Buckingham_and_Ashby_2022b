% This code draws heatmaps showing the rate of onset of resistance for
% different levels of pathogen mortality and sterility virulence. It also
% plots contours which show the proportion of the host population which is
% resistant.

% Set up parameters:
a0=5;
c1a=0.1;
c2a=-2;
b0=0.1;
c1b=0.1;
c2b=-2;
h=0.75;
beta=0.5;
delta=0.25;
q=1;
gamma=0;

% Vary mortality virulence (alpha) and sterility virulence (1-f)
alphavec_vals=linspace(0,1,51);
oneminusfvec_vals=linspace(0,0.99,50);
all_combinations=combvec(alphavec_vals,oneminusfvec_vals);
alphavec=all_combinations(1,:);
oneminusfvec=all_combinations(2,:);
varvec=alphavec;
ymax=1;

% Set up vectors to use later:
CSS1=NaN(length(varvec),1);
resprop1=NaN(length(varvec),1);
CSS2=NaN(length(varvec),1);
resprop2=NaN(length(varvec),1);
CSS3=NaN(length(varvec),1);
resprop3=NaN(length(varvec),1);
CSS4=NaN(length(varvec),1);
resprop4=NaN(length(varvec),1);
CSS5=NaN(length(varvec),1);
resprop5=NaN(length(varvec),1);
CSS6=NaN(length(varvec),1);
resprop6=NaN(length(varvec),1);

% For each combination of virulence parameters, we need to find the CSS
% value of the rate of onset of resistance and the corresponding resistant
% proportion.
parfor i=1:length(varvec)
    disp(i)
    
    alpha=alphavec(i);
    f=1-oneminusfvec(i);
    version=1;
    [CSSvec1,othervec1,~,respropvec1]=find_singstrats_function_heatmap(a0,c1a,c2a,b0,c1b,c2b,f,h,beta,alpha,delta,q,gamma,version,ymax);
    version=2;
    [CSSvec2,othervec2,~,respropvec2]=find_singstrats_function_heatmap(a0,c1a,c2a,b0,c1b,c2b,f,h,beta,alpha,delta,q,gamma,version,ymax);
    version=3;
    [CSSvec3,othervec3,~,respropvec3]=find_singstrats_function_heatmap(a0,c1a,c2a,b0,c1b,c2b,f,h,beta,alpha,delta,q,gamma,version,ymax);
    version=4;
    [CSSvec4,othervec4,~,respropvec4]=find_singstrats_function_heatmap(a0,c1a,c2a,b0,c1b,c2b,f,h,beta,alpha,delta,q,gamma,version,ymax);
    version=5;
    [CSSvec5,othervec5,~,respropvec5]=find_singstrats_function_heatmap(a0,c1a,c2a,b0,c1b,c2b,f,h,beta,alpha,delta,q,gamma,version,ymax);
    version=6;
    [CSSvec6,othervec6,~,respropvec6]=find_singstrats_function_heatmap(a0,c1a,c2a,b0,c1b,c2b,f,h,beta,alpha,delta,q,gamma,version,ymax);

    % We intend to plot heatmaps like these only in cases where the only
    % singular strategy is a single CSS. In other cases, we get these
    % warnings:
    if othervec1==1
        disp("Not just a CSS when version=1 and i=" +i)
    end
    if othervec2==1
        disp("Not just a CSS when version=2 and i=" +i)
    end
    if othervec3==1
        disp("Not just a CSS when version=3 and i=" +i)
    end
    if othervec4==1
        disp("Not just a CSS when version=4 and i=" +i)
    end
    if othervec5==1
        disp("Not just a CSS when version=5 and i=" +i)
    end
    if othervec6==1
        disp("Not just a CSS when version=6 and i=" +i)
    end
    
    % Create vectors of CSS values and resistant proportions:
    CSS1(i)=CSSvec1;
    resprop1(i)=respropvec1;
    CSS2(i)=CSSvec2;
    resprop2(i)=respropvec2;
    CSS3(i)=CSSvec3;
    resprop3(i)=respropvec3;
    CSS4(i)=CSSvec4;
    resprop4(i)=respropvec4;
    CSS5(i)=CSSvec5;
    resprop5(i)=respropvec5;
    CSS6(i)=CSSvec6;
    resprop6(i)=respropvec6;

end

% Change this data into matrices which can be plotted as heatmaps and
% contours:
matrix1=NaN(length(alphavec_vals),length(oneminusfvec_vals));
respropmat1=NaN(length(alphavec_vals),length(oneminusfvec_vals));
for j=1:length(oneminusfvec_vals)
    for i=1:length(alphavec_vals)
        matrix1(i,j)=CSS1(length(alphavec_vals)*(j-1)+i);
        respropmat1(i,j)=resprop1(length(alphavec_vals)*(j-1)+i);
    end
end
matrix2=NaN(length(alphavec_vals),length(oneminusfvec_vals));
respropmat2=NaN(length(alphavec_vals),length(oneminusfvec_vals));
for j=1:length(oneminusfvec_vals)
    for i=1:length(alphavec_vals)
        matrix2(i,j)=CSS2(length(alphavec_vals)*(j-1)+i);
        respropmat2(i,j)=resprop2(length(alphavec_vals)*(j-1)+i);
    end
end
matrix3=NaN(length(alphavec_vals),length(oneminusfvec_vals));
respropmat3=NaN(length(alphavec_vals),length(oneminusfvec_vals));
for j=1:length(oneminusfvec_vals)
    for i=1:length(alphavec_vals)
        matrix3(i,j)=CSS3(length(alphavec_vals)*(j-1)+i);
        respropmat3(i,j)=resprop3(length(alphavec_vals)*(j-1)+i);
    end
end
matrix4=NaN(length(alphavec_vals),length(oneminusfvec_vals));
respropmat4=NaN(length(alphavec_vals),length(oneminusfvec_vals));
for j=1:length(oneminusfvec_vals)
    for i=1:length(alphavec_vals)
        matrix4(i,j)=CSS4(length(alphavec_vals)*(j-1)+i);
        respropmat4(i,j)=resprop4(length(alphavec_vals)*(j-1)+i);
    end
end
matrix5=NaN(length(alphavec_vals),length(oneminusfvec_vals));
respropmat5=NaN(length(alphavec_vals),length(oneminusfvec_vals));
for j=1:length(oneminusfvec_vals)
    for i=1:length(alphavec_vals)
        matrix5(i,j)=CSS5(length(alphavec_vals)*(j-1)+i);
        respropmat5(i,j)=resprop5(length(alphavec_vals)*(j-1)+i);
    end
end
matrix6=NaN(length(alphavec_vals),length(oneminusfvec_vals)); 
respropmat6=NaN(length(alphavec_vals),length(oneminusfvec_vals));
for j=1:length(oneminusfvec_vals)
    for i=1:length(alphavec_vals)
        matrix6(i,j)=CSS6(length(alphavec_vals)*(j-1)+i);
        respropmat6(i,j)=resprop6(length(alphavec_vals)*(j-1)+i);
    end
end

% Make the plots:
subplot(2,3,1)
imagesc(oneminusfvec_vals,alphavec_vals,matrix1)
ylabel('Mortality virulence, $\alpha$','interpreter','latex')
title('Constant fecundity cost when resistant','interpreter','latex','fontsize',12)
set(gca,'YDir','normal')
colorbar
hold on
contour(oneminusfvec_vals,alphavec_vals,respropmat1,[0.1,0.4,0.5,0.55,0.6,0.65,0.7],'ShowText','on')
colormap(jet)
caxis([0,0.5])
text(0,0.96,'A','fontsize',20)

subplot(2,3,4)
imagesc(oneminusfvec_vals,alphavec_vals,matrix2)
title('Constant mortality cost when resistant','interpreter','latex','fontsize',12)
set(gca,'YDir','normal')
colorbar
hold on
contour(oneminusfvec_vals,alphavec_vals,respropmat2,[0.1,0.4,0.5,0.55,0.6,0.65,0.7],'ShowText','on')
colormap(jet)
caxis([0,0.5])
text(0,0.96,'D','fontsize',20)

subplot(2,3,2)
imagesc(oneminusfvec_vals,alphavec_vals,matrix3)
title('Whole-life reproduction trade-off','interpreter','latex','fontsize',12)
set(gca,'YDir','normal')
colorbar
hold on
contour(oneminusfvec_vals,alphavec_vals,respropmat3,[0.6,0.65,0.7,0.75],'ShowText','on')
colormap(jet)
caxis([0,0.5])
text(0,0.96,'B','fontsize',20)

subplot(2,3,5)
imagesc(oneminusfvec_vals,alphavec_vals,matrix4)
xlabel('Sterility virulence, $1-f$','interpreter','latex')
title('Whole-life mortality trade-off','interpreter','latex','fontsize',12)
set(gca,'YDir','normal')
colorbar
hold on
contour(oneminusfvec_vals,alphavec_vals,respropmat4,[0.6,0.65,0.7,0.75],'ShowText','on')
colormap(jet)
caxis([0,0.5])
text(0,0.96,'E','fontsize',20)

subplot(2,3,3)
imAlpha=ones(size(matrix5));
imAlpha(isnan(matrix5))=0;
imagesc(oneminusfvec_vals,alphavec_vals,matrix5,'AlphaData',imAlpha)
title('Pre-resistance reproduction trade-off','interpreter','latex','fontsize',12)
set(gca,'YDir','normal')
colorbar
hold on
contour(oneminusfvec_vals,alphavec_vals,respropmat5,[0.6,0.65,0.7,0.75],'ShowText','on')
colormap(jet)
caxis([0,0.5])
text(0,0.96,'C','fontsize',20)

subplot(2,3,6)
imAlpha=ones(size(matrix6));
imAlpha(isnan(matrix6))=0;
imagesc(oneminusfvec_vals,alphavec_vals,matrix6,'AlphaData',imAlpha)
title('Pre-resistance mortality trade-off','interpreter','latex','fontsize',12)
set(gca,'YDir','normal')
colorbar
ylabel('Rate of onset of resistance, $\zeta$','interpreter','latex')
hold on
contour(oneminusfvec_vals,alphavec_vals,respropmat6,[0.6,0.65,0.7,0.75],'ShowText','on')
colormap(jet)
caxis([0,0.5])
text(0,0.96,'F','fontsize',20)
