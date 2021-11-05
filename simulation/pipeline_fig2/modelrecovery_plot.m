function  modelrecovery_plot(CONFcorr,CONFacc,PerceptMod)

% redblue colormap
m = 64;
m1 = m*0.5;
r = (0:m1-1)'/max(m1-1,1);
g = r;
r = [r; ones(m1,1)];
g = [g; flipud(g)];
b = flipud(r);
 
cmap = [r g b];

% Belief recovery
f=figure(1);hold on
set(f,'name','P (simulated model | fit model) - belief trajectory')
for m = 1:length(PerceptMod)
    subplot(1,3,m);
    FM  = mean(CONFcorr.(PerceptMod{m}),3);
    IM  = FM./repmat(sum(FM),3,1);% inversion matrix
    IM(isnan(IM)) = 0; % 0 divisions causes nan
    IM = round(IM*1000)/1000;
    t   = imageTextMatrix(IM,{'State','Item','Comb'},{'State','Item','Comb'},[-1 1]);
    colormap(cmap)
    
    hold on;
    [l1, l2] = addFacetLines(zeros(3));
    set(t, 'fontsize', 30,'linewidth',3)
    set(gca,'TickLength',[0 0],'fontsize', 30)
    
    title(PerceptMod{m},'fontsize', 20,'fontweight', 'bold')
end
f.Position = [521 663 1800 500];

% Accuracy recovery
f=figure(2);hold on
set(f,'name','P (simulated model | fit model) - model trajectory')
for m = 1:length(PerceptMod)
    subplot(1,3,m);
    FM  = mean(CONFacc.(PerceptMod{m}),3);
    IM  = FM./repmat(sum(FM),3,1);% inversion matrix
    IM(isnan(IM)) = 0;% 0 divisions causes nan
    IM = round(IM*1000)/1000;
    t   = imageTextMatrix(IM,{'State','Item','Comb'},{'State','Item','Comb'},[-1 1]);
    colormap(cmap)
    
    hold on;
    [l1, l2] = addFacetLines(zeros(3));
    set(t, 'fontsize', 30,'linewidth',3)
    set(gca,'TickLength',[0 0],'fontsize', 30)
    
    title(PerceptMod{m},'fontsize', 20,'fontweight', 'bold')
end
f.Position = [528 89 1818 500];