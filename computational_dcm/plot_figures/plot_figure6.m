% function to reproduce figure 6
% -------------------------------
datatable   = readtable('minimum_dataset.csv');

% what to run
group_name  = {'NE','PTSD-','PTSD+'};
reg_name    = 'wHIP';


% plot info
grpcol      = {};
grpcol{1} = [153 204 255]./255;
grpcol{2} = [255 153 255]./255;
grpcol{3} = [204 0 102]./255;


figure(1); hold on

for gI = 1:3
    
    gid        = find(ismember(datatable.Group,group_name{gI}));
    X          = [datatable.(sprintf('Predictive_%s',reg_name)),datatable.(sprintf('Reactive_%s',reg_name))];% hippocampus
    X          = X(gid,:);
    
    imbalance = []; % in degree
    for s = 1:size(X,1)
        [F_R, theta,F_Ry,F_Rx] = forcesresultant(X(s,1),X(s,2));
        imbalance(s,1) = theta;
    end
    
    h = plotSpread(imbalance,'distributionColors',grpcol{gI},'distributionMarkers','o','xvalue',gI);
    M = rad2deg(circ_mean(deg2rad(imbalance)));% compute mean angle using circular toolbox function
    
    set(h{1},'Markersize',15,'markerfacecolor',grpcol{gI},'markeredgecolor',[0 0 0],'LineWidth',3)
    plot(gI-0.5,M,'o','Markersize',20,'markerfacecolor',[1 1 1],'markeredgecolor',grpcol{gI},'LineWidth',6)
    
    
    % bootstrap CI
    bootd1 = [];
    nboot = 1000;
    for b = 1:nboot
        
        % disp
        fprintf('performing bootstrap #%d\n',b)
        
        % boot sample
        d1 = imbalance;
        n1 = size(d1,1);
        r1 = randsample(n1,n1,1);
        
        % mean dir
        meandir1 = rad2deg(circ_mean(deg2rad(d1(r1,:))));
        
        % store
        bootd1 = [bootd1;meandir1];
        
    end
    
    ci      = [];
    alp     = .05;
    pct1    = 100*alp;
    pct2    = 100-pct1;
    
    % d1
    lower   = prctile(bootd1,pct1,1);
    upper   = prctile(bootd1,pct2,1);
    ci      = [lower;upper];
    line([gI-0.5 gI-0.5],[ci(1) ci(2)],'LineWidth',3,'color',[0 0 0])
end
set(gca,'ytick',[-180 -90 0 90 180])
set(gca,'xtick',[])
