datatable = readtable('minimum_dataset.csv');

% psycho-pathological dimension
Label   = {'Intrusion','Avoidance','Anhedonia','AnxArousal','Mood','DysArousal','Depression','SAnxiety'};
dim     = {};
dim{1}  = [1]; % Intrusion
dim{2}  = [2]; % Avoidance
dim{3}  = [4 6 8]; % Anxiety
dim{4}  = [3 5 7]; % Affect

% what to run
group_name  = {'PTSD-','PTSD+'};
reg_name    = 'wHIP';
typecor     = 'spearman';
alpha       = .05;
tail        = 'left';

% plot info
sLabel      = {'Intrusion','Avoidance','Anxiety','Affect'};
col         = [116,169,207; ... % Intrusion
    5,112,176; ... % Avoidance
    250,159,181; ... % Anxiety
    153,216,201;]./255;% Mood

nd          = length(sLabel);
figure; hold on;
shft = [0 0.2];ls={'--','-'};
subplot(1,2,1);hold on
for gI = 1:2

    gid     = find(ismember(datatable.Group,group_name{gI}));
    % loop trough dimension
    C = [];
    for i = 1:length(dim)

        X       =  datatable.(sprintf('Predictive_%s',reg_name))-datatable.(sprintf('Reactive_%s',reg_name));% hippocampus
        Y       = [];

        for ii = 1:length(dim{i})
            Y = [Y,datatable.(Label{dim{i}(ii)})];
        end

        x={};y={};
        x{1} = X(gid,:);y{1} = Y(gid,:);
        [CI_within, CI_between, es_within,es_between, p_within,p_between,creal] = bootcorr_figure5(x,y,typecor,alpha,tail,0);
        C(i) = creal;

        % bar
        ci = CI_within;
        line([i+shft(gI) i+shft(gI)],[ci(1) ci(2)],'LineWidth',2,'color',[0 0 0])
        plot(i+shft(gI),creal,'o','markersize',25,'markerfacecolor',col(i,:),'markeredgecolor',[0 0 0],'LineWidth',2);

    end
    %line
    plot([1:nd]+shft(gI),C,ls{gI},'LineWidth',2,'color',[0 0 0])

end
set(gca,'xtick',[1:nd],'xticklabel',sLabel)
ylim([-0.8 0.6])
xlim([0 nd+1])


% scatter plot avoidance
subplot(1,2,2);hold on;

X       =  datatable.(sprintf('Predictive_%s',reg_name))-datatable.(sprintf('Reactive_%s',reg_name));% hippocampus
Y       = [];
i       = 2;
for ii = 1:length(dim{i})
    Y = [Y,datatable.(Label{dim{i}(ii)})];
end
gid     = find(ismember(datatable.Group,group_name{2}));

X = X(gid,:);
Y = Y(gid,:);

P               = [X Y];
ot              = bivariate_outliers(zscore(P));% fonction from robust-correlation toolbox
P(ot>0,:)       = nan;
P(isnan(sum(P,2)),:) = [];

for p = 1:size(P,1)
    plot(P(p,1),P(p,2),'o','markersize',10,'markerfacecolor',col(1,:),'markeredgecolor','none');
end
[p,S]   = polyfit(P(:,1),P(:,2),1);
x1      = linspace(min(P(:,1)), max(P(:,1)));
y1      = polyval(p,x1);
plot(x1,y1,'LineWidth',3,'color',[0 0 0])
ylim([-0.1 0.3])
