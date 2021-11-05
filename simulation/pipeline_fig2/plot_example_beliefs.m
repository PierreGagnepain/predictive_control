function plot_example_beliefs(example_data)

f = figure;hold on;
plot(example_data.state,'color',[13,17,18]./255,'LineWidth',1.5);hold on;
plot(example_data.item,'color',[86,105,105]./255,'LineWidth',1.5);hold on;
plot(example_data.comb,'color',[24,118,158]./255,'LineWidth',1.5);hold on;

for t = 1:length(example_data.intrusions(:,1))
    plot(t,example_data.intrusions(t),'.','color',rand(1,3),'MarkerSize',10);hold on;
end

xlabel('Trials')
yticks([0 1])
yticklabels([{'Intrusion'};{'Non-intrusion'}])
ytickangle(90)
title('Trajectories of beliefs')
f.Position = [0 0 900 400];
