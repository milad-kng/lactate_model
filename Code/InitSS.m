function params = InitSS(params)

set(0,'defaultAxesFontSize',14,'defaultAxesFontName','Times',...
    'defaultAxesFontSmoothing','on',...
    'defaultLineLineWidth',1.3,'defaultAxesLineWidth',1.3,'defaultFigureColor','w')
if (params.glb.optSS == 1)
    % Volume Fractions
    params.vol.rea = 0.8;
    params.vol.rep = 0.444;
    params.vol.rca = 0.022;
    params.vol.VC = 0.0055;
    % Michaelis-Menten Parameters:
    params.MCT.EP.Vm = 1;% 3.1;%0.004-1%1;
    params.MCT.EP.Km = 2.2;%0.2-0.7-0.8
    % params.MCT.AE.Vm = 2;%0.1;%5;%37;%5;
    params.MCT.AE.Km = 28;%9.9-28-39
    params.MCT.AB.Vm = 0.5;%0.53%0.1
    params.MCT.AB.Km = 5.1;%1.9
    params.MCT.BA.Vm = 0.4;
    params.MCT.BA.Km = 5.1;%3.5-5.1-7.7%5.1
    % Blood Comp.
    params.Lac_j = 0.7;%0.88 0.8+-0.08 by Blood lactate is an important energy source for the human brain
    restt = [0.94,0.7,1,0.63,0.9,1.1,1.04,0.73];%1.5 excluded
    lacj = mean(restt);
    stdlacj = std(restt);
    params.CBF.F0 = 0.012;
    % LDH
    params.LDH.ast.Vn = 0.08;
    params.LDH.ast.Kn = 1;
    params.LDH.ka = 100;
    params.LDH.N1p = 18;
    % CBF
    %% Timing
    params.L = 30;%params.exp.p(end,1);
    params.dt = 0.005;
    params.CBF.t1 = 3;
    params.CBF.tend = 6;
    params.ft = 0:params.dt:params.L-params.dt;
    params.CBF.rep = 1;
%% Steady State: Optimization, visulization and selection %%
    numofsim = params.glb.numofsim;
    % vp+,kp+,vp-,kp-,lac_e,lac_c,va+,lac_p,lac_a
    bounds(1,:) = [0.5,0.03,0.24,0.3,0.2,0.66-2*0.11,0.5,0.8,0.8];%,0.05,0.05,0.5];
    bounds(2,:) = [60,0.07,28,8.5,5,0.66+2*0.11,70,1.2,1.2];%,0.3,0.3,1.5];
    % bounds(1,:) = [0,10,0.03*params.LDH.N1p/(8.5),0.3,0.2,0.66-2*0.11,...
    %                20,0.901-10*0.024,0.931-10*0.02,0.3];
    % bounds(2,:) = [10,28,0.07*params.LDH.N1p/(0.3),8.5,5,0.66+2*0.2,...
    %                50,0.901+10*0.024,0.931+10*0.02,1.5];
    params.dyn.mode = 0;
    params = glblopt(params,bounds,numofsim);
    save(strcat('params_Vm_AE_',num2str(params.MCT.AE.Vm),'Km_AE_',...
        num2str(params.MCT.AE.Km),'.mat'),'params')
else % Loading data
    if(params.glb.optSS==0)
        fig_type = params.glb.MCT4.fig_type;
        switch fig_type
            case "typical"
                a = load(strcat('params_new_nov2021_Vm_AE_',...
                    num2str(params.MCT.AE.Vm),'Km_AE_',...
                num2str(params.MCT.AE.Km),'.mat'),'params');
                params = a.params;
                params.dyn.mode = 0;
                params.glb.optSS = 0;
                params = SelectSS(params);
            case "MCT4"% Full stat figures
                load(strcat('params_VAE_total.mat'),'params');
                VAE_vector = params.glb.MCT4.VAE;
                for index = 1:length(VAE_vector)
                    params.MCT.AE.Vm = VAE_vector(index);
                    params.dyn.mode = 0;
                    params.glb.optSS = 0;
                    params.glb.MCT4.loop_num = index;
                    params.glb.MCT4.fig_type = fig_type;
                    params = SelectSS(params);
                    rates{index} = params.glb.MCT4.vals{1,index}.flows;
                end
                params.glb.MCT4.rates = rates;
                params.glb.MCT4.VAEs = VAE_vector;
                plotMeanFlow(params);
                params.glb.MCT4.fig_type = 'MCT4_summary';
                plotMeanFlow(params);
%                 param_dist_VAE(params,1);
            case 'only param selection'

        end
    end
end
end
%% -----------------------------------------------------------------------
function params = SelectSS(params)
params.glb.method = "ICA";
if strcmp(params.glb.method,"ICA")
    s = 1;
    F_values = params.glb.fval;
    else
if strcmp(params.glb.method,"FG")
    s = 2;
    F_values = sum(params.glb.fval,2);
end
end

sols = reshape(params.glb.sols(:,s,:),size(params.glb.sols(:,s,:),1),...
    size(params.glb.sols(:,s,:),3));
j = 0;
% F_scalar = sum(params.glb.fval,2);
for m = 1:size(sols,1)
    if (F_values(m,:)<10^(-25))
        j=j+1;
        fval_new(j) = params.glb.fval(m);
        sols_new(j,:) = sols(m,:);
        solica_ind(j,:) = m;
   end
end
ss_ica = [sols_new(:,5:6) sols_new(:,8:9)];
dss_mean = ss_ica - mean(ss_ica);
dssPA = abs(dss_mean(:,3))+abs(dss_mean(:,4));
[c, d] = min(dssPA);
[c, d] = min(fval_new);

[a, b] = min(abs(dss_mean)); 
disp(d)
params.glb.best_a = sols_new(b(4),:);
params.glb.best_p = sols_new(b(3),:); 
if (mod(length(sols_new(:,1)),2)==0)
    params.glb.bestSS = sols_new(length(sols_new(:,1))/2,:);%sols_new(d,:);%should be d
else
    params.glb.bestSS = sols_new((1+length(sols_new(:,1)))/2,:);
end
params.glb.bestSS = sols_new(end,:);
    % P:389, 445, 109
% A:4,195
% params.glb.ind = solica_ind;
params.pyr_P = params.glb.bestSS(8)/params.LDH.N1p;
params.pyr_A = params.glb.bestSS(9)/params.LDH.ka;
% fluxes = plotMeanFlow(sols_new,params);
% params.glb.flows = fluxes;
params.glb.solica = sols_new;
params.glb.solF = fval_new;
for j = 1:size(sols_new,1)
    flx(j,:) = syslac(sols_new(j,:),params,2);
end
flows = [flx(:,1:5) flx(:,6)-flx(:,7) flx(:,8)-flx(:,9)];
flows(:,5) = flows(:,5)*params.vol.VC;
params.glb.flows = flows;
end
%% 
function plotMeanFlow(params)
titlef = {'V_E_P';'V_A_E';'V_A_C-V_C_A';'V_C_A_P';...
          'Pyramidal cell';'Astrocyte'};
switch params.glb.MCT4.fig_type
    case "typical"
        flows = params.glb.flows;
        flowsm = mean(flows);
        flowst = std(flows);
        names = ["V_P^P^r^o^d","K_P^P^r^o^d","V_P^C^o^n^s",...
                 "K_P^C^o^n^s","Lac_E^S^S","Lac_C^S^S",...
                 "V_A^P^r^o^d","Lac_P^S^S","Lac_A^S^S"];
        figure
        for i = 1:9
            subplot(3,3,i)
            histogram(solica(:,i),20)
            title(names(i)) 
        end
        figure
        subplot(1,2,1)
        errorbar(1:4,[flowsm(1:2) flowsm(3)-flowsm(4) flowsm(5)],[flowst(1:2) (flowst(3)-flowst(4)) flowsm(5)],...
            '*','MarkerSize',15,'CapSize',10,'LineWidth',1.5)
        set(gca,'xtick',1:4,'xticklabel',titlef(1:4))
        ylabel(['net exchange rate (mM/min)',' (mean',char(177),'std)'])
        xlabel('exchanges')
        xlim([0 4])
        title('(a)')
        % grid on
        subplot(1,2,2)
        errorbar(1:2,flowsm(6:7),flowst(6:7),'*','MarkerSize',15,'CapSize',10,'LineWidth',1.5)
        set(gca,'xtick',1:2,'xticklabel',titlef(5:6))
        ylabel(['production-consumption difference (mM/min)',' (mean',char(177),'std)'])
        xlabel('compartment')
        xlim([0 3])
        title('(b)')
        % grid on
        
        % boxplots
        figure
        subplot(1,3,1)
        boxplot([flows(:,1:2)],1:2)
        set(gca,'TickLabelInterpreter', 'tex');
        set(gca,'xtick',1:2,'xticklabel',titlef(1:2))
        ylabel('net exchange rate (mM/min)')
        xlabel('exchanges')
        xlim([0 3])
        title('(a)')
        % grid on

        subplot(1,3,2)
        boxplot(flows(:,6:7),1:2)
        set(gca,'TickLabelInterpreter', 'tex');
        set(gca,'xtick',1:2,'xticklabel',titlef(5:6))
        ylabel('production-consumption difference (mM/min)')
        xlabel('compartment')
        xlim([0 3])
        title('(b)')
        % grid on

        subplot(1,3,3)
        boxplot([flows(:,3)-flows(:,4) flows(:,5)],1:2)
        set(gca,'TickLabelInterpreter', 'tex');
        set(gca,'xtick',1:2,'xticklabel',titlef(3:4))
        ylabel('net exchange rate (mM/min)')
        xlabel('exchanges')
        xlim([0 3])
        title('(c)')
        % grid on
        % resting concentrations comparison
        figure
        errorbar(1:2,[mean(solica(:,8)) mean(solica(:,9))],[std(solica(:,8)) std(solica(:,9))],...
            '*','MarkerSize',15,'CapSize',10,'LineWidth',1.5)
        set(gca,'xtick',1:2,'xticklabel',titlef(5:6))
        ylabel(['resting concentrations (mM)', '(mean',char(177),'std)'])
        xlabel('compartment')
        xlim([0 3])
        title('(b)')
        % grid on
        [h ,p] = ttest2(solica(:,8),solica(:,9));
        if (p<0.05)
            disp('Null hypothesis rejected')
        end
    case "MCT4"
        rates = params.glb.MCT4.rates;
        VAE_vector = params.glb.MCT4.VAE;
        n = length(VAE_vector);
        groups = [];
        rates_combined = [];
        for i = 1:n
            C(i,:) = blanks(4);
            C(i,1:length(num2str(VAE_vector(i),2))) =...
             num2str(VAE_vector(i),2);
            if (i == 2 || i==3)
                C(2,1:length(num2str(VAE_vector(2),2))) = '0.5 ';
                C(3,1:length(num2str(VAE_vector(3),2))) = '1   ';
            end
            temp = repmat(C(i,:),length(rates{1,i}),1);
            flows_median(i,:) = median(rates{1,i});
            flows_mean(i,:) = mean(rates{1,i});
            groups = [groups;temp];
            rates_combined = [rates_combined;rates{1,i}];
            num_positive(i,:) = sum(rates{1,i}>0);
            num_negative(i,:) = sum(rates{1,i}<0);
        end
        ratioPN = num_positive./num_negative;
        f1 = figure('Name','Rates','Position', [10 10 1200 600]);
        Names = ["V_E_N","V_A_C","J_A"];
        for ii=1:2:5
%             ax(ii).YAxis.TickLabelFormat = '%.4e';
%             f(ii) = figure('Name','Ratio and median');
%             ax_curve = gca;
%             ax_curve.YLabel.String = Names (ii);
%             ax_curve.XLabel.String = 'V_A_E^m';
%             ax_curve.YAxis.Exponent        = 0;
%             ax_curve.YAxis.TickLabelFormat = '%.4e';
            if (ii==1)
                flows = rates_combined(:,ii);
                median_vec = flows_median (:,ii);
                mean_vec = flows_mean(:,ii);
            else
                if(ii==3)
                   flows = rates_combined(:,3)-rates_combined(:,4);
                   median_vec = flows_median (:,3) - flows_median (:,4);
                   mean_vec = flows_mean(:,3) - flows_mean(:,4);
                end    
                if(ii==5)
                   flows = rates_combined(:,7);
                   median_vec = flows_median (:,7);
                   mean_vec = flows_mean (:,7);
                end
            end
            figure(f1)
            for sub_idx = ii:ii+1
                ax(sub_idx) = subplot(3,2,sub_idx);
                ax(sub_idx).YLabel.String = Names(floor(ii/2)+1);
                ax(sub_idx).XLabel.String = 'V_A_E^m';
                ax(sub_idx).FontWeight = 'bold';
                ax(sub_idx).YAxis.Exponent = 0;

                switch mod(sub_idx,2)
                    case 1
                        hold(ax(sub_idx),'on')
                        ax(sub_idx).YLabel.String = Names(floor(ii/2)+1);
                        ax(sub_idx).XAxis.TickLabelFormat = '%.2g';
                        boxplot(ax(sub_idx),flows,groups)%,'plotstyle','compact'
                        set(ax(sub_idx),'TickLabelInterpreter', 'tex','ylim'...
                            ,[min(flows)*0.9 max(flows)*1.05]);
                        ax(sub_idx).XTickLabelRotation = 90;
                        pos_2 = ax(sub_idx).Position(2);
                        pos_3 = ax(sub_idx).Position(3);
                        pos_4 = ax(sub_idx).Position(4);
                        disp(ax(sub_idx).Position)
                    case 0
%                         ax(sub_idx).XTick = VAE_vector; % I should change xtick labels for 
                                % second column to strings
                                % of this
                                % vector
                        set(ax(sub_idx),'TickLabelInterpreter', 'tex')
                        hold(ax(sub_idx),'on')
%                         ax(sub_idx).XTickLabelRotation = 90;
                        ax(sub_idx).XAxis.TickLabelFormat = '%.2g';
                        plot(ax(sub_idx),VAE_vector,mean_vec,'Color','k')
%                         hold(ax(sub_idx),'on')
%                         plot(ax(sub_idx),VAE_vector,median_vec)
                        box(ax(sub_idx),'on')
                        set(ax(sub_idx),'TickLabelInterpreter', 'tex','ylim'...
                            ,[0.9*min([min(mean_vec),min(median_vec)]) 
                        1.05*max([max(mean_vec),max(median_vec)])]);
                        ax(sub_idx).Position(4) = pos_4;
                        ax(sub_idx).Position(2) = pos_2;
                        ax(sub_idx).Position(3) = pos_3/2;
                        disp(ax(sub_idx).Position)
%                         legend('mean')
                        if (sub_idx==2)
                            hold(ax(sub_idx),'on')
                            plot(ax(sub_idx),VAE_vector,zeros(length(VAE_vector)))
                        end
                end
            end
        end
        
    case 'MCT4_summary'
        rates = params.glb.MCT4.rates;
        n = 20;
        VAE_vector = params.glb.MCT4.VAEs;
        for i = 1:n
            C(i,:) = blanks(4);
            C(i,1:length(num2str(VAE_vector(i),2))) =...
             num2str(VAE_vector(i),2);
            if (i == 2 || i==3)
                C(2,1:length(num2str(VAE_vector(2),2))) = '0.5 ';
                C(3,1:length(num2str(VAE_vector(3),2))) = '1   ';
            end
            temp = repmat(C(i,:),length(rates{1,i}),1);
%            flows_median(i,:) = median(rates{1,i});
            flows_mean(i,:) = mean(rates{1,i});
        end
        f1 = figure('Name','Rates','Position', [10 10 1300 300]);
        Names = ["V_E_P (mM/min)","V_A_C (mM/min)","J_A (mM/min)"];
        for ii=1:2:5
            if (ii==1)
%                 median_vec = flows_median (:,ii+1);
                mean_vec = flows_mean(:,ii);
            else
                if(ii==3)
%                    median_vec = flows_median (:,3) - flows_median (:,4);
                   mean_vec = flows_mean(:,3) - flows_mean(:,4);
                end    
                if(ii==5)
%                    median_vec = flows_median (:,7);
                   mean_vec = flows_mean (:,7);
                end
            end
            figure(f1)
            sub_idx = (ii+1)/2;
                ax(sub_idx) = subplot(1,3,sub_idx);
                ax(sub_idx).YAxis.Exponent = 0;
                set(ax(sub_idx),'TickLabelInterpreter', 'tex')
                hold(ax(sub_idx),'on')
%                         ax(sub_idx).XTickLabelRotation = 90;
                ax(sub_idx).XAxis.TickLabelFormat = '%.2g';
                plot(ax(sub_idx),VAE_vector,mean_vec,'color','k')
                if (sub_idx==1)
                    hold(ax(sub_idx),'on')
                    plot(ax(sub_idx),VAE_vector,zeros(length(VAE_vector)))
                end
%                 hold(ax(sub_idx),'on')
%                 plot(ax(sub_idx),VAE_vector,median_vec)
                box(ax(sub_idx),'on')
                set(ax(sub_idx),'TickLabelInterpreter', 'tex','xlim'...
                    ,[0.1 8])
                if(sub_idx==2)
                    xlabel('V_A_E^m','fontweight','bold');
                end
                ylabel (Names(sub_idx),'fontweight','bold'); 
                ax(sub_idx).XAxis.TickLabelFormat = '%.2g';
                end
            end
end



%%

