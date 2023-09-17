function plot_borgperf_pcycneg(rxnList, Stat, set_idx, outfp, neg_idx)
%Function to plot perfomance measures on a list of true reaction
%assignments, and statistics from borgifier
rng('default')
P_idx=set_idx(rxnList(set_idx)==Stat.bestMatchIndex(set_idx));
N_idx=set_idx(rxnList(set_idx)~=Stat.bestMatchIndex(set_idx));

%fill up  negatives with to create equal set
Tot_N=find(rxnList==-1 & ismember(Stat.bestMatchIndex, neg_idx));
Add_N_idx=randsample(Tot_N, length(P_idx)-length(N_idx));
N_idx=[N_idx;Add_N_idx];
set_idx=[set_idx; Add_N_idx];
%plot count density histograms with equal bins
close all
edges = linspace(0, 1, 11);
hist1 = histogram(Stat.bestMatch(P_idx),edges, 'FaceColor', 'blue', 'EdgeColor', 'none');
hold on; % Hold the current plot

hist2 = histogram(Stat.bestMatch(N_idx),edges, 'FaceColor', 'red', 'EdgeColor', 'none');

% Add labels and a legend
xlabel('Score');
ylabel('Density');
legend('P', 'N');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 7 5]);
if ~isfolder(fileparts(outfp))
    mkdir(fileparts(outfp))
end
saveas(gcf, [outfp, 'pcycneg_dens.svg'])
close gcf


% plot rc curve
%sort reaction List by scores
[~, sidx]=sort(Stat.bestMatch(set_idx),"descend");
set_idx=set_idx(sidx);
TPR=[];
FPR=[];
precision=[];
thres=[];
for i=1:20
    TPR=[TPR; sum(ismember(set_idx(1:ceil((length(set_idx)/20)*i)),P_idx))/length(P_idx)];
    FPR=[FPR; sum(ismember(set_idx(1:ceil((length(set_idx)/20)*i)),N_idx))/length(N_idx)];
    precision=[precision; sum(ismember(set_idx(1:ceil((length(set_idx)/20)*i)),P_idx))/ceil((length(set_idx)/20)*i)];
    thres=[thres; Stat.bestMatch(set_idx(ceil((length(set_idx)/20)*i)))];
end
tiledlayout(1,2)
nexttile
plot(FPR,TPR, '-', 'LineWidth', 2);
title('ROC curve')

nexttile
plot(thres, precision, '-', 'LineWidth', 2)
xlabel('Threshhold');
ylabel('Precicision');
xlim([0,1])
ylim([0,1])
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 10 5]);
saveas(gcf, [outfp, 'pcycneg_meas.svg'])
close gcf
