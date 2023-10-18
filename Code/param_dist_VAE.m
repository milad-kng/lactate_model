function []= param_dist_VAE(params,p)
set(0,'defaultAxesFontSize',14,'defaultAxesFontName','Times',...
    'defaultAxesFontSmoothing','on',...
    'defaultLineLineWidth',1.3,'defaultAxesLineWidth',1.3,'defaultFigureColor','w',...
    'defaultAxesFontWeight','Bold')
n = length(params.glb.MCT4.VAE);
P_pos = [];
P_neg = [];
k = 500;
P_neg_random = zeros(k,9);
P_pos_random = zeros(k,9);
params.dyn.artery.VAE = [];
for i = 1:n
    VAE = params.glb.MCT4.VAE(i);
    vals = params.glb.MCT4.vals{1,i};
    flows = vals.flows(:,p);
    flows_cap = vals.flows(:,3) - vals.flows(:,4);
    sum(flows_cap<0)
    P_neg_loop = vals.solica(flows<0,:);
    P_pos_loop = vals.solica(flows>0,:);
    P_neg = [P_neg;vals.solica(flows<0,:)];
    P_pos = [P_pos;vals.solica(flows>0,:)];
    rand_idx_neg = randsample(size(vals.solica(flows<0,:),1),k/n);
%     fprintf('neg = %d\n',rand_idx_neg)
    rand_idx_pos = randsample(size(vals.solica(flows>0,:),1),k/n);
%     fprintf('pos = %d\n',rand_idx_pos)
    P_neg_random((i-1)*(k/n)+1:i*(k/n),:) = P_neg_loop(rand_idx_neg,:);
    P_pos_random((i-1)*(k/n)+1:i*(k/n),:) = P_pos_loop(rand_idx_pos,:);
    params.dyn.artery.VAE = [params.dyn.artery.VAE; repmat(VAE,k/n,1)];
end

%% PLot resting
P_neg_params = [P_neg(:,1:4) P_neg(:,7)];
P_pos_params = [P_pos(:,1:4) P_pos(:,7)];
title_str = ["V_N^P^r^o^d","K_N^P^r^o^d","V_N^C^o^n^s",...
                 "K_N^C^o^n^s","V_A^P^r^o^d"];
figure('Position', [30 30 900 150])
for j = 1:5
        subplot(1,5,j)
%         figure('Position', [30 30 200 200])
        errorbar(1,mean(P_neg_params(:,j)),std(P_neg_params(:,j)),'.'...
            ,'MarkerSize',15,'CapSize',10,'LineWidth',1.5);
        hold on;
        errorbar(2,mean(P_pos_params(:,j)),std(P_pos_params(:,j)),'.'...
            ,'MarkerSize',15,'CapSize',10,'LineWidth',1.5);
        set(gca,'TickLabelInterpreter', 'tex');
        set(gca,'xtick',1:2,'xticklabel','')
        ylim(gca, 'padded')
        
        y1=ylim(gca);
        yticks([y1(1) 0.5*(y1(1)+y1(2)) y1(2)])
        yticklabels({num2str(y1(1)) num2str(0.5*(y1(1)+y1(2))) ...
                    num2str(y1(2))});
        title(title_str(j),'fontweight','bold');
%         [h, p_val] = ttest2(P_neg_params(:,j),P_pos_params(:,j));
%         ylabel(num2str(p_val))
        xlim([0 3])
        ytickformat('%.2f');
%         ax=gca;
%         outer = ax.OuterPosition;
%         ax.OuterPosition = outer + [0.05 0 -0.05 0];
        if (j==5)
            legend('Lactate release','Lactate uptake')
        end
end
conc_neg = [P_neg(:,8) P_neg(:,9) P_neg(:,5) P_neg(:,6)];
conc_pos = [P_pos(:,8) P_pos(:,9) P_pos(:,5) P_pos(:,6)];

% title('Resting concentrations for negative regime')
% ---------------Negative Regime------------------ %
figure('Position', [0 0 320 120])%450
subplot(1,2,1)
%neuron
errorbar(1,mean(conc_neg(:,1)),std(conc_neg(:,1)),...
    '.','MarkerSize',15,'CapSize',10,'LineWidth',1.5,'Color','r');
xlabel('Resting conentrations for lactate release behavior','fontweight','bold')
hold on;
%astrocyte
errorbar(2,mean(conc_neg(:,2)),std(conc_neg(:,2)),...
    '.','MarkerSize',15,'CapSize',10,'LineWidth',1.5,'Color',winter(1));
%extracel
errorbar(3,mean(conc_neg(:,3)),std(conc_neg(:,3)),...
    '.','MarkerSize',15,'CapSize',10,'LineWidth',1.5);
hold on;
% capillary
errorbar(4,mean(conc_pos(:,4)),std(conc_pos(:,4)),...
    '.','MarkerSize',15,'CapSize',10,'LineWidth',1.5);
% legend('neuron','astrocyte','extracellular','capillary','Location','northeastoutside')

% legend('neuron','astrocyte')
set(gca,'xtick',1:2,'xticklabel','')
ylim(gca, 'padded')
y1=ylim(gca);
yticks([y1(1) 0.5*(y1(1)+y1(2)) y1(2)])
yticklabels({num2str(y1(1)) num2str(0.5*(y1(1)+y1(2))) ...
                    num2str(y1(2))});
xlim([0 5])
ytickformat('%.2f');
% ax=gca;
% outer = ax.OuterPosition;
% ax.OuterPosition = outer + [0.05 0 -0.05 0];
% ----------------- Positive regime ---------------------%
subplot(1,2,2)
%figure('Position', [30 30 200 200])
% neuron
errorbar(1,mean(conc_pos(:,1)),std(conc_pos(:,1)),...
    '.','MarkerSize',15,'CapSize',10,'LineWidth',1.5,'Color','r');
xlabel('Resting conentrations for lactate release behavior','fontweight','bold')
hold on;
%astrocyte
errorbar(2,mean(conc_pos(:,2)),std(conc_pos(:,2)),...
    '.','MarkerSize',15,'CapSize',10,'LineWidth',1.5,'Color',winter(1));
%extracel
errorbar(3,mean(conc_pos(:,3)),std(conc_pos(:,3)),...
    '.','MarkerSize',15,'CapSize',10,'LineWidth',1.5);
%xlabel('Lactate uptake regime','fontweight','bold')
hold on;
% capillary
errorbar(4,mean(conc_pos(:,4)),std(conc_pos(:,4)),...
    '.','MarkerSize',15,'CapSize',10,'LineWidth',1.5);
legend('neuron','astrocyte','extracellular','capillary','Location','northeastoutside')
set(gca,'xtick',1:4,'xticklabel','')
ylim(gca, 'padded')
y1=ylim(gca);
yticks([y1(1) 0.5*(y1(1)+y1(2)) y1(2)])
yticklabels({num2str(y1(1)) num2str(0.5*(y1(1)+y1(2))) ...
                    num2str(y1(2))});
xlim([0 5])
ytickformat('%.2f');
% ax=gca;
% outer = ax.OuterPosition;
% ax.OuterPosition = outer + [0.05 0 -0.05 0];
%% Effect of artery lactate
if (params.dyn.artery.mode == 1)
    lac_art = 5*params.Lac_j;
    params.dyn.artery.neg.samples = P_neg_random;
    params.dyn.artery.pos.samples = P_pos_random;
    params.dyn.artery.lac_art = lac_art;
    % Simulate negative samples
    params = state_computation(params,P_neg_random);
    % plotting_aux(params)
    diff_min_neg = params.dyn.artery.diff_min;
    diff_max_neg = params.dyn.artery.diff_max;
    states_neg = params.dyn.artery.states;
    interm_states_neg = params.dyn.artery.interm_states;
    % prepare to save
    params.dyn.artery.neg.states = states_neg;
    params.dyn.artery.neg.interm_states = interm_states_neg;
    params.dyn.artery.neg.diff_min = diff_min_neg;
    params.dyn.artery.neg.diff_max = diff_max_neg;
    % simulate positive values
    params = state_computation(params,P_pos_random);
    % plotting_aux(params)
    diff_min_pos = params.dyn.artery.diff_min;
    diff_max_pos = params.dyn.artery.diff_max;
    states_pos = params.dyn.artery.states;
    interm_states_pos = params.dyn.artery.interm_states;
    % prepare to save
    params.dyn.artery.pos.states = states_pos;
    params.dyn.artery.pos.interm_states = interm_states_pos;
    params.dyn.artery.pos.diff_min = diff_min_pos;
    params.dyn.artery.pos.diff_max = diff_max_pos;
    save('params_final_paper_revised','params')
else
    load('params_final_paper_revised','params')
    interm_states_pos = params.dyn.artery.pos.interm_states;
    diff_min_pos = params.dyn.artery.pos.diff_min;
    diff_max_pos = params.dyn.artery.pos.diff_max;
    interm_states_neg = params.dyn.artery.neg.interm_states;
    diff_min_neg = params.dyn.artery.neg.diff_min;
    diff_max_neg = params.dyn.artery.neg.diff_max;
    % --------------- Visualization ----------%
    %-------VAE----------%
    figure
    th = 2000;
    interm_states_neg_clean = interm_states_neg;%(1:2,interm_states_neg(1,:,1:2)<th,:);
    interm_states_pos_clean = interm_states_pos;%(1:2,interm_states_pos(1,:,1:)<th,:);
%     diff_min_pos = diff_min_pos(1:2,interm_states_neg(1,:,1)<th,1);
%     diff_min_neg = diff_min_neg(1:2,interm_states_neg(1,:,1)<th,1);
    [h, p_val] = ttest2(interm_states_pos_clean(1,:,2),interm_states_pos_clean(2,:,2));
    fprintf('\nV_A_E p-value for positive regime = %d\n',p_val)
    subplot(1,2,1)
    boxplot((squeeze(interm_states_pos_clean(1:2,:,2)))',1:2,...
        'Notch','off','ExtremeMode','compress','Symbol','','Labels',...
        {'Lac_j_0','10 Lac_j_0'});
    h1 = findobj(gcf, 'tag', 'Upper Adjacent Value');
    set(gca,'TickLabelInterpreter', 'tex');
    ylim([-0.005 0.09])
    ax = gca;
    ax.YAxis.Exponent = -2;
%     ylim([-1 1.1*h1(1,:).YData(1)])
%     title('Resting uptake')
%     title('V_A_E changes as a function of artery lactate increase in positive regime')
%     ylabel('V_A_E / V_A_E^S^S')
    ylabel('max (V_A_E) - V_A_E^S^S')
    [h, p_val] = ttest2(interm_states_neg_clean(1,:,2),interm_states_neg_clean(2,:,2));
    fprintf('\nV_A_E p-value for negative regime  = %d \n',p_val)
    subplot(1,2,2)
    boxplot((squeeze(interm_states_neg_clean(1:2,:,2)))',1:2,...
        'Notch','off','ExtremeMode','compress','Symbol','','Labels',...
        {'Lac_j_0','10 Lac_j_0'});
    h1 = findobj(gcf, 'tag', 'Upper Adjacent Value');
    h2 = findobj(gcf, 'tag', 'Lower Whisker');
    set(gca,'TickLabelInterpreter', 'tex');
    ylim([-0.005 0.09])
    ax = gca;
    ax.YAxis.Exponent = -2
%     ylim([-1 1.1*h1(1,:).YData(1)])
%     title('Resting release')
%     title('V_A_E changes as a function of artery lactate increase in negative regime')
%     ylabel('V_A_E / V_A_E^S^S')
    ylabel('max (V_A_E) - V_A_E^S^S')
    saveas(gcf,'C:\Users\milad\OneDrive - Concordia University - Canada\Research_new\Research\Lactate\Paper\imgs\box_AE.emf')

    %-------VEN Pos and negative mean-------

    figure
    [h, p_val] = ttest2(interm_states_pos_clean(1,:,1),interm_states_pos_clean(2,:,1));
    fprintf('\nV_E_N p-value for positive regime  = %d \n',p_val)
    subplot(1,2,1)
    boxplot((squeeze(interm_states_pos_clean(1:2,:,1)))',1:2,...
        'Notch','off','ExtremeMode','compress','Symbol','','Labels',...
        {'Lac_j_0','10 Lac_j_0'});
    h1 = findobj(gcf, 'tag', 'Upper Adjacent Value');
    h2 = findobj(gcf, 'tag', 'Lower Adjacent Value');
    set(gca,'TickLabelInterpreter', 'tex');
    ylim([0.00 0.09])
    ax = gca;
    ax.YAxis.Exponent = -2;
%     ylim([-1 1.1*h1(1,:).YData(1)])
%     title('Resting uptake')
%     title('V_E_N changes as a function of artery lactate increase in positive regime')
%     ylabel('V_E_N / V_E_N^S^S')
    ylabel('max (V_E_N) - V_E_N^S^S')
    [h, p_val] = ttest(interm_states_neg_clean(1,:,1),interm_states_neg_clean(2,:,1));
    fprintf('\nV_E_N p-value for negative regime  = %d \n',p_val)
    subplot(1,2,2)
    boxplot((squeeze(interm_states_neg_clean(1:2,:,1)))',1:2,...
        'Notch','off','ExtremeMode','compress','Symbol','','Labels',...
        {'Lac_j_0','10 Lac_j_0'});
    h1 = findobj(gcf, 'tag', 'Upper Adjacent Value');
    h2 = findobj(gcf, 'tag', 'Lower Adjacent Value');
    set(gca,'TickLabelInterpreter', 'tex');
    ylim([0.00 0.09])
    ax = gca;
    ax.YAxis.Exponent = -2;
%     ylim([-1 1.1*h1(1,:).YData(1)])
%     title('Resting release')
%     legend('positive resting state','negative resting state')
%     title('V_E_N changes as a function of artery lactate increase in negative regime')
%     ylabel('V_E_N / V_E_N^S^S')
    ylabel('max (V_E_N) - V_E_N^S^S')
    saveas(gcf,'C:\Users\milad\OneDrive - Concordia University - Canada\Research_new\Research\Lactate\Paper\imgs\box_EN.emf')
    %-------------J_N----------------%
    figure
    [h, p_val] = ttest2(interm_states_pos_clean(1,:,6),interm_states_pos_clean(2,:,6));
    fprintf('\nJ_N p-value for positive regime  = %d \n',p_val)
    subplot(1,2,1)
    boxplot((squeeze(interm_states_pos_clean(1:2,:,6)))',1:2,...
        'Notch','off','ExtremeMode','compress','Symbol','','Labels',...
        {'Lac_j_0','10 Lac_j_0'});
    h1 = findobj(gcf, 'tag', 'Upper Adjacent Value');
    h2 = findobj(gcf, 'tag', 'Lower Adjacent Value');
    set(gca,'TickLabelInterpreter', 'tex');
    ylim([-0.03 0.0])
    ax = gca;
    ax.YAxis.Exponent = -2;

%     ytickangle(90)
%     ylim([-1 1.1*h1(1,:).YData(1)])
%     title('Resting uptake')
%     title('V_E_N changes as a function of artery lactate increase in positive regime')
%     ylabel('V_E_N / V_E_N^S^S')
    ylabel('max (J_N) - J_N^S^S')
    [h, p_val] = ttest(interm_states_neg_clean(1,:,6),interm_states_neg_clean(2,:,6));
    fprintf('\nJ_N p-value for negative regime  = %d \n',p_val)
    subplot(1,2,2)
    boxplot((squeeze(interm_states_neg_clean(1:2,:,6)))',1:2,...
        'Notch','off','ExtremeMode','compress','Symbol','','Labels',...
        {'Lac_j_0','10 Lac_j_0'});
    h1 = findobj(gcf, 'tag', 'Upper Adjacent Value');
    h2 = findobj(gcf, 'tag', 'Lower Adjacent Value');
    set(gca,'TickLabelInterpreter', 'tex');
    ylim([-0.03 0])
    ax = gca;
    ax.YAxis.Exponent = -2;
%     ylim([-1 1.1*h1(1,:).YData(1)])
%     title('Resting release')
%     legend('positive resting state','negative resting state')
%     title('V_E_N changes as a function of artery lactate increase in negative regime')
%     ylabel('V_E_N / V_E_N^S^S')
    ylabel('max (J_N) - J_N^S^S')
    saveas(gcf,'C:\Users\milad\OneDrive - Concordia University - Canada\Research_new\Research\Lactate\Paper\imgs\box_JN.emf')
%-------------J_A----------------%
    figure
    [h, p_val] = ttest2(interm_states_pos_clean(1,:,7),interm_states_pos_clean(2,:,7));
    fprintf('\nJ_N p-value for positive regime  = %d \n',p_val)
    subplot(1,2,1)
    boxplot((squeeze(interm_states_pos_clean(1:2,:,7)))',1:2,...
        'Notch','off','ExtremeMode','compress','Symbol','','Labels',...
        {'Lac_j_0','10 Lac_j_0'});
    h1 = findobj(gcf, 'tag', 'Upper Adjacent Value');
    h2 = findobj(gcf, 'tag', 'Lower Adjacent Value');
    set(gca,'TickLabelInterpreter', 'tex');
    ylim([-0.00005 0.0008])
%     ylim([-1 1.1*h1(1,:).YData(1)])
%     set(gca, 'yticklabel', get(gca, 'ytick'));
%     ax = gca;
%     ax.YAxis.Exponent = -2;
%     title('Resting uptake')
%     title('V_E_N changes as a function of artery lactate increase in positive regime')
%     ylabel('V_E_N / V_E_N^S^S')
    ylabel('max (J_A) - J_A^S^S')
    [h, p_val] = ttest(interm_states_neg_clean(1,:,7),interm_states_neg_clean(2,:,7));
    fprintf('\nJ_N p-value for negative regime  = %d \n',p_val)
    subplot(1,2,2)
%     set(gca, 'yticklabel', get(gca, 'ytick'));
    boxplot((squeeze(interm_states_neg_clean(1:2,:,7)))',1:2,...
        'Notch','off','ExtremeMode','compress','Symbol','','Labels',...
        {'Lac_j_0','10 Lac_j_0'});
    h1 = findobj(gcf, 'tag', 'Upper Adjacent Value');
    h2 = findobj(gcf, 'tag', 'Lower Adjacent Value');
    set(gca,'TickLabelInterpreter', 'tex');
 ylim([-0.00005 0.0008])
%     ylim([-1 1.1*h1(1,:).YData(1)])
%     title('Resting release')
%     legend('positive resting state','negative resting state')
%     title('V_E_N changes as a function of artery lactate increase in negative regime')
%     ylabel('V_E_N / V_E_N^S^S')
    ylabel('max (J_A) - J_A^S^S')
%     set(gca, 'yticklabel', get(gca, 'ytick'));
    saveas(gcf,'C:\Users\milad\OneDrive - Concordia University - Canada\Research_new\Research\Lactate\Paper\imgs\box_JA.emf')
    
    % --------------Undershoot----------
    figure
    subplot(1,2,1)
    [h, p_val] = ttest2(diff_min_pos(1,:,1),diff_min_pos(2,:,1));
    fprintf('\nundershoot p-value for positive regime  = %d \n',p_val)
    boxplot((squeeze(diff_min_pos(1:2,:,1)))',1:2,...
        'Notch','off','ExtremeMode','compress','Symbol','','Labels',...
        {'Lac_j_0','10 Lac_j_0'});
    h1 = findobj(gcf, 'tag', 'Upper Adjacent Value');
    h2 = findobj(gcf, 'tag', 'Lower Adjacent Value');
    set(gca,'TickLabelInterpreter', 'tex');
    ylim([-0.2 -0.02])
    ax = gca;
    ax.YAxis.Exponent = -2;
%     ylim([0.9*h2(2,:).YData(1) 1.2*h1(1,:).YData(1)])
%     title('Resting uptake')
%     legend('positive resting state','negative resting state')
%     title('V_E_N changes as a function of artery lactate increase in negative regime')
    ylabel('min (Lac_N) - Lac_N^S^S')

    subplot(1,2,2)
    [h, p_val] = ttest2(diff_min_neg(1,:,1),diff_min_neg(2,:,1));
    fprintf('\nundershoot p-value for negative regime  = %d \n',p_val)
    boxplot((squeeze(diff_min_neg(1:2,:,1)))',1:2,...
        'Notch','off','ExtremeMode','compress','Symbol','','Labels',...
        {'Lac_j_0','10 Lac_j_0'});
    h1 = findobj(gcf, 'tag', 'Upper Adjacent Value');
    h2 = findobj(gcf, 'tag', 'Lower Adjacent Value');
    set(gca,'TickLabelInterpreter', 'tex');
    ylim([-0.2 -0.02])
    ax = gca;
    ax.YAxis.Exponent = -2;
%     ylim([0.9*h2(2,:).YData(1) h1(1,:).YData(1)])
%     title('Resting release')
%     legend('positive resting state','negative resting state')
%     title('V_E_N changes as a function of artery lactate increase in negative regime')
    ylabel('min (Lac_N) - Lac_N^S^S')
    saveas(gcf,'C:\Users\milad\OneDrive - Concordia University - Canada\Research_new\Research\Lactate\Paper\imgs\undershoot.emf')

%     h1 = findobj(gcf, 'tag', 'Outliers')
    %-----------Overshoot--------%
    figure
    subplot(1,2,1)
    [h, p_val] = ttest(diff_max_pos(1,:,1),diff_max_pos(2,:,1));
    fprintf('\novershoot p-value for positive regime  = %d \n',p_val)
    boxplot((squeeze(diff_max_pos(1:2,:,1)))',1:2,...
        'Notch','off','ExtremeMode','compress','Symbol','','Labels',...
        {'Lac_j_0','10 Lac_j_0'});
    h1 = findobj(gcf, 'tag', 'Upper Adjacent Value');
    h2 = findobj(gcf, 'tag', 'Lower Adjacent Value');
    set(gca,'TickLabelInterpreter', 'tex');
    ylim([-0.002 0.02])
    ax = gca;
    ax.YAxis.Exponent = -2;
%     ylim([0.995*h2(2,:).YData(1) 1.005*h1(1,:).YData(1)])
%     ylim([-0.001 1.1*h1(1,:).YData(1)])

%     title('Resting uptake')
%     legend('positive resting state','negative resting state')
%     title('V_E_N changes as a function of artery lactate increase in negative regime')
    ylabel('max (Lac_N) - Lac_N^S^S')
%     set(gca, 'yticklabel', get(gca, 'ytick'));
    
    [h, p_val] = ttest(diff_max_neg(1,:,1),diff_max_neg(2,:,1));
    fprintf('\novershoot p-value for negative regime  = %d \n',p_val)
    subplot(1,2,2)
    boxplot((squeeze(diff_max_neg(1:2,:,1)))',1:2,...
        'Notch','off','ExtremeMode','compress','Symbol','','Labels',...
        {'Lac_j_0','10 Lac_j_0'});
    h1 = findobj(gcf, 'tag', 'Upper Adjacent Value');
    h2 = findobj(gcf, 'tag', 'Lower Adjacent Value');
    set(gca,'TickLabelInterpreter', 'tex');
    ylim([-0.002 0.02])
    ax = gca;
    ax.YAxis.Exponent = -2;
%     ylim([-0.001 1.1*h1(1,:).YData(1)])
%     title('Resting release')
%     legend('positive resting state','negative resting state')
%     title('V_E_N changes as a function of artery lactate increase in negative regime')
    ylabel('max (Lac_N) - Lac_N^S^S')
%     set(gca, 'yticklabel', get(gca, 'ytick'));
    saveas(gcf,'C:\Users\milad\OneDrive - Concordia University - Canada\Research_new\Research\Lactate\Paper\imgs\overshoot.emf')
%     ytickformat('%.4f');
end
%%%%%%%%%-------------------Computing states-----------------%%%%%%%%%%%%%
function params = state_computation(params,samples)
f5 = figure('Name','inputs','Position', [10 10 1100 180]);
tau = [0.4 0.4 0.4 0.4 1 1 0.4 0.4];
params.dyn.mode = 2;
params.dyn.func = ["SHS","SHS","DE","SHS"];
params.L = 30;%params.exp.p(end,1);
params.dt = 0.005;
params.CBF.t1 = 3;
params.CBF.tend = 6;
params.ft = 0:params.dt:params.L-params.dt;
params.CBF.rep = 1;
% states = zeros(length(samples),length(params.ft),4);
% % interm_states = zeros(2,length(samples),length(params.ft),7);
% diff_min = zeros(2,length(samples),4);
% diff_max = zeros(2,length(samples),4);
for i = 1:length(samples)
    param_set = samples(i,:);
    ss(i,:) = [param_set(8),param_set(5),param_set(9),param_set(6)];
    params.glb.bestSS = param_set;
    params.pyr_P = ss(i,1)/params.LDH.N1p;
    params.pyr_A = ss(i,3)/params.LDH.ka;
    v = [params.pyr_P*0.8,params.pyr_A,params.CBF.F0*0.6,...
         params.Lac_j*10];
    params.MCT.AE.Vm = params.dyn.artery.VAE(i);
    [t,X1] = ode15s(@(t,x)sys4D(t,x,params,0,[tau v],f5),params.ft,ss(i,:));
    % [V_ep;V_ae;V_ac;V_ca;V_cap;Pp;Cp;Pa;Ca];
    flows1 = sys4D(t',X1',params,1,[tau v],f5);
    v = [params.pyr_P*0.8,params.pyr_A,params.CBF.F0*0.6,...
        params.Lac_j];
    [t,X0] = ode15s(@(t,x)sys4D(t,x,params,0,[tau v],f5),params.ft,ss(i,:));
    flows0 = sys4D(t',X0',params,1,[tau v],f5);
%     states(i,:,:) = (X1-X0)./ss;
    states(i,:,:) = (X1-X0)./ss(i,:);
    ss_flows = (flows1(:,5));
    if(ss_flows(1)<0)
        L = -1;
    else
        L = 1;
    end
    max_flows = max(flows0(:,params.ft>3 & params.ft<6),[],2);
%     interm_states(1,i,:) = abs((max_flows)./((ss_flows)));
    interm_states(1,i,:) = (max_flows)-(ss_flows);
    max_flows = max(flows1(:,params.ft>3 & params.ft<6),[],2);
%     interm_states(2,i,:) = abs((max_flows)./((ss_flows)));
    interm_states(2,i,:) = (max_flows)-(ss_flows);

%     diff_min(1,i,:) = (min(X0,[],1))./ss(i,:);
%     diff_min(2,i,:) = (min(X1,[],1))./ss(i,:);
%     diff_max(1,i,:) = (max(X0,[],1))./ss(i,:);
%     diff_max(2,i,:) = (max(X1,[],1))./ss(i,:);

    diff_min(1,i,:) = min(X0,[],1) - ss(i,:);
    diff_min(2,i,:) = min(X1,[],1) - ss(i,:);
    diff_max(1,i,:) = max(X0,[],1) - ss(i,:);
    diff_max(2,i,:) = max(X1,[],1) - ss(i,:);
    fprintf('%d\n',min(X1,[],1)-min(X0,[],1))
%     figure;plot(t,X1(:,1))
%     hold on;plot(t,X0(:,1))
%     figure;plot(t,X1(:,2))
%     hold on;plot(t,X0(:,2))
    fprintf('\n next:\n')
end
params.dyn.artery.states = states;
params.dyn.artery.interm_states = interm_states;
params.dyn.artery.t = t;
params.dyn.artery.diff_min = diff_min;
params.dyn.artery.diff_max = diff_max;
params.dyn.artery.ss = ss;
%% %%%%%%%-------------------Plotting paper results-----------------%%%%%%%%%%%%%
function [] = plotting_aux(params)

f1 = figure('Name','desired concentrations and rates');
figure (f1)
Names = ["Lac_P (mM)";"Lac_E (mM)";"V_E_P (mM/min)";"V_A_E (mM/min)"];
for a=1:4
   ax(a) = subplot(2,2,a);
   ax(a).YLabel.String = Names(a);
   ax(a).XLabel.String = 'time (min)';
end
t1 = params.CBF.t1;
tend = params.CBF.tend;
dt = 3;
states = params.dyn.artery.states;
states_interm = params.dyn.artery.interm_states;
t = params.dyn.artery.t;
% colors = colormap(winter(100));
for i = 1:size(params.dyn.artery.states,1)
     variable_to_plot = [reshape(states(i,:,1),length(params.ft),1),reshape(states(i,:,2),length(params.ft),1),...
        reshape(states_interm(i,:,1),length(params.ft),1),reshape(states_interm(i,:,2),length(params.ft),1)];
    for r = 1:4
        figure(f1)
        hold (ax(r),'on')
        box(ax(r),'on')
        colormap(ax(r),spring)
        plot(ax(r),t,variable_to_plot(:,r));
        if(i == length(params.dyn.artery.states))
            ylim(ax(r), 'auto')%'tight'
            g_x = t1:dt:tend; 
            g_y = get(ax(r),'ylim');
            for l=1:length(g_x)
                plot(ax(r),[g_x(l) g_x(l)],g_y,'k:') %y grid lines  
                hold on
                plot(ax(r),[t1 tend],[g_y(2) g_y(2)],'k:')
            end
            xticks(ax(r),[t1:dt:tend (tend+2*dt:dt*2:params.L)]);
            xticklabels(ax(r),['0','3',string(tend+dt:dt*2:params.L)]);
        end
    end
end
% for kk=2:2
%     hold(ax(kk),'on')
%     legend(legnds) 
% end





