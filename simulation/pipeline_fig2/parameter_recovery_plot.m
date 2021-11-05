function [simfit_corr] = parameter_recovery_plot(simulated, fitted,PerceptMod,ObsMod,parameters,nsim,real_sim_diff, real_sim_corr)


idsup = 1:nsim;

for m = 1:length(PerceptMod)
    for par = 1:length(parameters{m})
        for ob = 1:length(ObsMod)
            
            simu = []; 
            fitt = [];
            
            simu = squeeze(simulated.(PerceptMod{m}).(ObsMod{ob}).(parameters{m}{par}));
            fitt = squeeze(fitted.(PerceptMod{m}).(ObsMod{ob}).(parameters{m}{par}));
            
            if m == 1 && ob == 1
                fitt(fitt==-3)=nan;                
            end
            
            for i = 1:length(idsup)
                
                df = real_sim_diff.(PerceptMod{m}).(ObsMod{ob})(i,:)';
                co = real_sim_corr.(PerceptMod{m}).(ObsMod{ob})(i,:)';
                
                idk = 1:100;%intersect(find(abs(df)<0.1),find(co>0.5));
                
                try
                    r = corr(simu(idsup(i),idk)',fitt(idsup(i),idk)', 'type','Spearman','rows','complete');
                catch
                    r = nan;
                end
                
              
                simfit_corr.(PerceptMod{m}).(ObsMod{ob}).(parameters{m}{par})(i,1) = r;
                
            end
        end
    end
end


%% PLOT PARAMETER RECOVERY

% Reshape simfit_corr to plot it
fscorr = [];

fscorr = [{simfit_corr.HGF.state.om},{simfit_corr.HGF.item.om},...
    {simfit_corr.RW.state.al},{simfit_corr.RW.item.al},...
    {simfit_corr.KF.state.expom},{simfit_corr.KF.item.expom},...
    {simfit_corr.KF.state.pi_u},{simfit_corr.KF.item.pi_u}];
    
% Plot specifications
col =[204,255,230
    204,255,230
    232,255,48
    232,255,48
    255,222,194
    255,222,194
    255,222,194
    255,222,194]./255;

xshift = 1:8;
sW      = 0.5;
sMM     = 0;
cn      = 0;
grpcol  = [153 204 255;255 153 255;204 0 102]./255;

fm = @(x)nanmean(x);

f = figure;hold on;
set(f,'Name','Parameter recovery')
for pp = 1:8
    d = [fscorr{pp}];

    xV      = [0.5+xshift(pp)];
    
    hold on;
    plotSpread(d,'xValues', xV,'spreadWidth',sW,...
        'showMM',sMM,'distributionColors',col(pp,:));
    hold on;
        % ci
    ci      = [];   
    ci(:,1) = bootci(2000,{fm, d(:,1)},'type','per');
 
   % line([0.8+xshift(pp) 0.8+xshift(pp)],[ci(1,1) ci(2,1)],'LineWidth',2,'color',[0 0 0])
   % plot([0.8+xshift(pp)],mean(d),'or','color',[0 0 0],'LineWidth',3)
   plot([0.8+xshift(pp)],nanmean(d),'.','MarkerEdgeColor',[0 0 0],'LineWidth',3)
    

end
set(findall(gca,'type','line'),'markerSize',18)
set(findall(gca,'type','dot'),'markerSize',36)
set(findall(gca,'type','dot'),'markerEdgeSize',2)
set(findall(gca,'type','dot'),'markerEdgeColor',[29 29 27])


set(gca,'xtick',1+[1:8],'xticklabel',{'Omega (HGF)',[],'Alpha (RW)',[],'Expom (KF)',[],'PiU (KF)',[]})
ylabel('Correlation')

f.Position = [100 100 900 500];

%% Report mean and standerd deviation
% 
% meanPR = []; stdPR = [];
% for e = 1:length(toeval)
%     for ob = 1:2
%         meanPR(e,ob) = mean(c_simfit{e}(:,ob));
%         stdPR(e,ob)  = std(c_simfit{e}(:,ob));
%     end
% end

