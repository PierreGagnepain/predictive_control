
function [Falsification, SimIntrusionPattern,real_sim_diff, real_sim_corr] = model_falsification_plot(sim_intrusions,PerceptMod,ObsMod,nsim,nsampling,nsess)

% real intrusion data across TNT sessions (average of whole population)
real_intrusions = [0.3886    0.3109    0.2656    0.2474];


for m = 1:length(PerceptMod)
    
    % response across TNT block
    intstate = [];intitem = [];intcomb = [];
    Bstate = [];Bitem = [];Bcomb = [];
    for s = 1:nsim
        
        for ob = 1:length(ObsMod)
            sim        = [];
            sim        = sim_intrusions{s}.(PerceptMod{m}).(ObsMod{ob});
            
            % store mean session intrusions across sampled parameters
            mean_sim.(PerceptMod{m}).(ObsMod{ob})(s,:) = mean(sim);
            
            for p = 1:nsampling
                
                x       = real_intrusions;
                y       = sim(p,:);
                
                % mean difference between real intrusions and simulated intrusions
                real_sim_diff.(PerceptMod{m}).(ObsMod{ob})(s,p)       = mean(x-y);
                
                % correlation between real intrusions and simulated intrusions
                real_sim_corr.(PerceptMod{m}).(ObsMod{ob})(s,p)       = corr(x',y');
            end
            SimIntrusionPattern.(PerceptMod{m}).(ObsMod{ob})(s,:) = mean(sim);
        end
    end
end

% Compute mean error across all the simulations
mean_diff = []; mean_corr = [];std_diff = [];std_corr = [];
for m = 1:length(PerceptMod)
    for ob = 1:length(ObsMod)
        mean_diff(m,ob) = nanmean(nanmean(real_sim_diff.(PerceptMod{m}).(ObsMod{ob})));
        mean_corr(m,ob) = nanmean(nanmean(real_sim_corr.(PerceptMod{m}).(ObsMod{ob})));
        
        std_diff(m,ob) = std(nanmean(real_sim_diff.(PerceptMod{m}).(ObsMod{ob})));
        std_corr(m,ob) = std(nanmean(real_sim_corr.(PerceptMod{m}).(ObsMod{ob})));
        
    end
end

analyses_name = {'Mean difference';'Mean difference';'Mean difference';...
    'Correlation';'Correlation';'Correlation'};

Falsification = table(analyses_name,repmat(PerceptMod',2,1),[mean_diff(:,1);mean_corr(:,1)],...
    [mean_diff(:,2);mean_corr(:,2)],[mean_diff(:,3);mean_corr(:,3)]);

Falsification.Properties.VariableNames = {'Analysis','Percept',ObsMod{1},ObsMod{2},ObsMod{3}};


% Plot real and simulated intrusion patterns
f = @(x)mean(x);

col =[13,17,18
    86,105,105
    24,118,158]./255;
spaceinc = [0 0.2 0.4];

P.markertype = {'s','v','p'}; % state,item,combined, marker plots

fp = figure; hold on;

for m = 1:length(PerceptMod)
    subplot(1,3,m); hold on;
    
    % Real intrusion pattern
    plot(real_intrusions, 'LineWidt',5,'color',[0.7 0 0]);hold on;
    plot(1:4,real_intrusions,'o','MarkerSize',15,'MarkerEdgeColor',[0.7 0 0],'LineWidth',3,'MarkerFaceColor',[1 1 1])
    
    for ob = 1:length(ObsMod)
        
        tp = mean(mean_sim.(PerceptMod{m}).(ObsMod{ob}));
        
        plot([1:4]+spaceinc(ob),tp, 'LineWidth',2.5,'color',col(ob,:));hold on;
        
        for ss = 1:nsess
            
            xc = ss+spaceinc(ob);
            
            D = mean_sim.(PerceptMod{m}).(ObsMod{ob})(:,ss);
            ci = std(D);%bootci(2000,f,mean_sim.(PerceptMod{m}).(ObsMod{ob})(:,ss));
            %line([xc xc],[mean(D)+ci mean(D)-ci],'Color',[0 0 0],'LineWidt',1.5)
            line([xc xc],[mean(D)+ci mean(D)-ci],'Color',col(ob,:),'LineWidt',2)
            
            plot(xc,tp(ss),P.markertype{ob},'MarkerSize',10,'MarkerEdgeColor',col(ob,:),'LineWidth',3,'MarkerFaceColor',col(ob,:))
        end
        
        % Simulated intrusion pattern
        % plot(ss,tp(ss),'o','MarkerSize',10,'LineWidt',5,'MarkerEdgeColor',P.col(ob,:))
        
        ylim([0 0.5])
        xlim([-0.5 5])
        set(gca,'xtick',1:4)
        set(gca,'ytick',[0 0.5])
        xlabel('Session')
        ylabel('Intrusion rating')
        hYLabel = get(gca,'ylabel');
        set(hYLabel,'rotation',90,'VerticalAlignment','middle')
        
    end
    title(PerceptMod{m},'fontweight','bold','fontsize',15)
end


set(fp,'Name','Model falsification')

fp.Position = [100 100 900 500];





