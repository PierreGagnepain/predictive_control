% function to reproduce figure 4
% -------------------------------
datatable   = readtable('minimum_dataset.csv');

group_name  = {'NE','PTSD-','PTSD+'};
region_name = {'rHIP','cHIP','PC','wHIP'};
cond_name   = {'Reactive','Predictive'};
condcol     = repmat([216 39 23;5 105 140],3,1)./255;
fm          = @(x)mean(x);
for reg = 1:length(region_name)
    
    figure(reg);hold on;
    
    % to plot
    inc = 0;
    TP = {};
    M = [];
    for gI = 1:3 % group
        gid = find(ismember(datatable.Group,group_name{gI}));

        for c = 1:2 % condition
            inc = inc + 1;
            dat     = datatable.(sprintf('%s_%s',cond_name{c},region_name{reg}));
            dat     = dat(gid);
            TP{inc} = dat;
            M(inc)  = mean(dat);
        end
        
    end
    
    % plot
    h = plotSpread(TP,'distributionColors',condcol,'distributionMarkers','o','spreadWidth',0.5);
    ch = get(h{3},'children');
    
    
    idh = [6:-1:1];
    for g = 1:6
        
        
        set(ch(idh(g)),'Markersize',10,'markerfacecolor',condcol(g,:),'markeredgecolor',[0 0 0],'LineWidth',3)
        
        % add line between mean
        meanpos = [g+0.3 g+1-0.3];
        if ~isEven(g)
            line(meanpos,[M(g) M(g+1)],'LineWidth',2,'color',[0 0 0])
            
            % add dot
            plot(meanpos(1),M(g),'o','Markersize',12,'markerfacecolor',[1 1 1],'markeredgecolor',condcol(g,:),'LineWidth',6)
            plot(meanpos(2),M(g+1),'o','Markersize',12,'markerfacecolor',[1 1 1],'markeredgecolor',condcol(g+1,:),'LineWidth',6)
            
            % add CI
            ci = bootci(2000,{fm, TP{g}},'type','per');
            line([meanpos(1) meanpos(1)],[ci(1) ci(2)],'LineWidth',2,'color',condcol(g,:))
            ci = bootci(2000,{fm, TP{g+1}},'type','per');
            line([meanpos(2) meanpos(2)],[ci(1) ci(2)],'LineWidth',2,'color',condcol(g+1,:))
            
        end
        
    end
    ylim([-4 4])
    set(gca,'xtick',[])
end