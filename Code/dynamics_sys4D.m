%% Dynamic
function [params] = dynamics_sys4D(params,multi)%purturb, loop is deleted
ss = [params.glb.bestSS(8),params.glb.bestSS(5),params.glb.bestSS(9),...
      params.glb.bestSS(6)];
init_val = [params.pyr_P,params.pyr_A,params.CBF.F0,params.Lac_j];
names = ["Pyr_N","Pyr_A","CBF","Lac_j"];
set(0,'defaultAxesFontSize',14,'defaultAxesFontName','Times',...
    'defaultAxesFontSmoothing','on',...
    'defaultLineLineWidth',1.3,'defaultAxesLineWidth',1.3,'defaultFigureColor','w')
if (multi == 0 || multi == 1 || multi==2)
    f1 = figure('Name','Concentrations','Position', [10 10 500 400]);
    f3 = figure('Name','rates','Position', [10 10 520 640]);
%     f4 = figure('Name','LDH activity');
    f5 = figure('Name','inputs','Position', [10 10 1100 180]);
    figure(f1)
    Names = ["Lac_N (mM)";"Lac_E (mM)";"Lac_A (mM)";"Lac_C (mM)"];
    for a=1:4
       ax(a) = subplot(2,2,a);
       ax(a).YLabel.String = Names(a);
       ax(a).YLabel.FontWeight = 'bold';
%        ax(a).XLabel.String = 'time (min)';
    end
    %set(ax,'FontWeight','bold','FontName','Calibri','FontSmoothing','on')
%     set(f1,'RendererMode','manual','DefaultFigureRenderer','opengl')
    figure(f3)
    Names = ["V_E_N (mM/min)","V_A_E (mM/min)","V_A_C (mM/min)",...
                "V_C_A_P (mM/min)","J_N (mM/min)","J_A (mM/min)"];
    for a=1:6
       axe(a) = subplot(3,2,a);
       axe(a).YLabel.String = Names(a);
       axe(a).YLabel.FontWeight = 'bold';
%        axe(a).XLabel.String = 'time (min)';
    end
figure(f5)    
coef = ones(1,4);
if(multi == 0)
    s = purturb;
    f2 = figure('Name','Input');
    figure(f2)
    ax2 = subplot(1,1,1);
for i = 1:length(loop)
    coef(s) = loop(i);
    v = coef.*init_val;
    [t,X]=NaiveRK4Sto(@(t,x)sys4D(t,x,v,params,0,func,f5),params.L,params.dt,0,ss);
    t = t(1,1:end-2);
    X = X(:,1:end-2);
    active = tdf(params,func(s),init_val(s),v(s));
    active = interp1(params.ft,active,t);
    flows = sys4D(t,X,v,params,1,func);
    if (max(active) == init_val(s))
        legnds = strcat(char(names(s)),'=',num2str(ceil(min(active)/init_val(s))),'*',char(names(s)),'0');
    else
        legnds = strcat(char(names(s)),'=',num2str(ceil(max(active)/init_val(s))),'*',char(names(s)),'0');
    end
    for j = 1:4
        figure(f1)
        hold (ax(j),'on')
        plot(ax(j),t,X(j,:),'DisplayName',legnds)
    end
    for k = 1:4
       figure(f3) 
       hold (axe(k),'on')
       if (k<=2 || k==4)
           plot(axe(k),t,flows(k,:),'DisplayName',legnds)
           legend (axe(k),'show')
       end
       if (k == 3)
           plot(axe(k),t,flows(k,:)-flows(k+1,:),'DisplayName',legnds)
           legend (axe(k),'show')
       end
    end
    figure(f4)
    for n = 1:2
        hold (ax4(n),'on')
        plot(ax4(n),t,flows(n+5,:),'DisplayName',legnds)
        legend (ax4(n),'show')
        grid (ax4(n),'on')
    end

    figure(f2)
    hold (ax2,'on')
    plot(t,active,'DisplayName',legnds)
end
figure(f1)
% suptitle(['Concentration (mM) vs. time(min)','type:',func(s)])
for j = 1:4
%     figure(f1)
%     hold (ax(j),'on')
%     plot(t,active*100,'--','DisplayName','Input curve')
    legend (ax(j),'show')

    h1 = line([params.CBF.t1 params.CBF.t1],[min(X(j,:)) max(X(j,:))]);
    h2 = line([params.CBF.tend params.CBF.tend],[min(X(j,:)) max(X(j,:))]);
    set([h1 h2],'Color','k','LineWidth',2)
% Add a patch
    patch([params.CBF.t1 params.CBF.tend params.CBF.tend params.CBF.t1],[min(X(j,:)) min(X(j,:)) max(X(j,:)) max(X(j,:))],'r')
    set([h1 h2],'Color','k','LineWidth',2)
end

figure(f2)
legend (ax2,'show')
ax2.Title.String = ['input:',names(s), func(s)];
params.dyn.X = X;
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(multi == 1)
    c = 5;
    m = 0;
    tau = params.dyn.tau;
    v = params.dyn.v;
%     for u = 1:4
%         if(strcmp(params.dyn.func(u),"")
%     end
%     v(1) = params.dyn.v(1)*params.pyr_P/(1-(exp((-(params.CBF.tend-params.CBF.t1))/tau(1))));
%     v(2) = params.dyn.v(2)*params.pyr_A/(1-(exp((-(params.CBF.tend-params.CBF.t1))/tau(3))));
    S = 4;
    for variable = linspace(params.Lac_j,params.Lac_j*10,c)
        m = m+1;
        v(S) = variable;
        [t,X] = ode15s(@(t,x)sys4D(t,x,params,0,[tau v],f5),params.ft,ss);
%         vec = X(:,1)/ss(1);
%         decay_timeP(m) = decay_calc(vec);
%         vec = X(:,3)/ss(3);
%         decay_timeA(m) = decay_calc(vec);
        flows = sys4D(t',X',params,1,[tau v],f5);
        init_val = [params.pyr_P,params.pyr_A,params.CBF.F0,params.Lac_j];
        init = init_val(S);
        features = [init,variable,tau(7:8)];
        form = params.dyn.func(S);
        active = tdf(params,features,form);
        legnds{m} = strcat(names(S),'=',num2str(max(active)/init_val(S)),names(S),'_0');
        plotting('flow');
%         plotting('LDH');
        plotting('conc');
    end
%     figure(f1)
%   suptitle(['Concentration (mM) vs. time(min)']);%,'type:',func(s)])
%     for j = 2:2
%         legend (ax(j),'show')
%     end
    params.dyn.X = X;
    params.dyn.time = t;
end
% -------------------------------------------------------
if(multi == 2)
    lacj = linspace(params.Lac_j*3,params.Lac_j*10,15);
    tau = params.dyn.tau;
    for i = 1:length(lacj)
        v = loop;
        v(3) = lacj(i);
        [t,X] = ode15s(@(t,x)sys4D(t,x,params,0,[tau v],f5),params.ft,ss);
        flows = sys4D(t,X,v,params,1,func);
        input =  tdf(params,'DEN',init_val(1),v(1));
        m2m = min(input)/params.pyr_P;
        init_val = [params.pyr_P,params.pyr_A,params.CBF.F0,params.Lac_j];
        legnds = strcat(char(names(1)),'=',num2str(m2m),char(names(1)),'_0');
        for s = 2:4
           legnds = strcat(legnds,newline,',',char(names(s)),'=',num2str(loop(s)/init_val(s)),char(names(s)),'_0');
        end
       legnds = strcat('Pyr_P=',num2str(m2m),'pyr_P_0');
        for k = 1:4
           figure(f3) 
           hold (axe(k),'on')
           if (k<=2 || k==4)
            plot(axe(k),t,flows(k,:),'DisplayName',legnds,'LineWidth',1.3)
            legend (axe(k),'show')
           else
               if (k == 3)
                   plot(axe(k),t,flows(k,:)-flows(k+1,:),'DisplayName',legnds,'LineWidth',1.3)
                   legend (axe(k),'show')
               end
           end
        end
        figure(f4)
        for n = 1:2
            hold (ax4(n),'on')
            plot(ax4(n),t,flows(n+5,:),'DisplayName',legnds,'LineWidth',1.3)
            legend (ax4(n),'show')
            grid(ax4(n),'on')
        end

        for j = 1:4
            figure(f1)
            hold (ax(j),'on')
            plot(ax(j),t,X(j,:),'DisplayName',legnds,'LineWidth',1.3)
            grid(ax(j),'on')
            legend(ax(j),'show')
        end
    end
    end
    % figure(f1)
    % suptitle(['concentration (mM) vs. time(min)']);%,'type:',func(s)])
    for j = 1:4
    %     figure(f1)
    %     hold (ax(j),'on')
    %     plot(t,active*100,'--','DisplayName','Input curve')
    %    legend (ax(j),'show')
    end
    end
%------------initial dip testing--------------%
if (multi == 3)% dynamics vs lacj
    c = 8;
    f5 = figure('Name','inputs');
%     ax=get(gca);
%     set(gca,'FontWeight','bold','FontSize',12,'FontName','Calibri','FontSmoothing','on')
%     ax_f = subplot(1,1,1);
    tau = params.dyn.tau;
    f6 = figure('Name','demand vs artery lac.');
    ax_j = subplot(1,1,1);
    lac_art = linspace(params.Lac_j,params.Lac_j*10,c);
    pyr_ast = linspace(params.pyr_A,params.pyr_A*5,c);
%     v_tend = v0 -v1*(1-(exp((-(params.CBF.tend-params.CBF.t1))/tau(1))));
    pyr_pyr = 0.98*params.pyr_P;%linspace(v_tend,v_tend,1) ;% 0.15
    for i = 1:length(pyr_pyr)
        for j = 1:length(pyr_ast)
%             f(j) = figure('Name',strcat('flows',pyr_ast(j)));
%             figure(f(j))
%             ax_f(1) = subplot(3,1,1);
%             ax_f(2) = subplot(3,1,2);
%             ax_f(3) = subplot(3,1,3);
            for k = 1:length(lac_art)
                v(1) = pyr_pyr(i);
                v(2) = pyr_ast(j);
                v(3) = params.CBF.F0*1.5;
                v(4) = lac_art(k);
                [t,X] = ode15s(@(t,x)sys4D(t,x,params,0,[tau v]),params.ft,ss);
                flows = sys4D(t',X',params,1,[tau v],f5);
                if ((ss(1) - min(X(:,1)))>=0)
                    demand(i,j,k) = ss(1) - min(X(:,1));
%                     demand_flows (i,j,k) = flows(1,1) - min() 
                else
                    demand(i,j,k) = 0;
                end
%                 if (k==1)
                input =  tdf(params,[init_val(4),v(3),tau(end-1:end)],'DE2');
                lgn = strcat('Lac_J=',num2str(max(input)/params.Lac_j),'Lac_J_0');
%                 hold(ax_f(1),'on')
%                 plot(ax_f(1),t,flows(1,:),'DisplayName',lgn)
%                 legend(ax_f(1),'show')
%                 ylabel(ax_f(1),'V_E_P')
%                 hold(ax_f(2),'on')
%                 plot(ax_f(2),t,flows(6,:)-flows(7,:))
% %                 legend(ax_f(2),'show')
%                 ylabel(ax_f(2),'P_P-C_P')
%                 hold(ax_f(3),'on')
%                 plot(ax_f(3),t,X(:,1))
% %                 legend(ax_f(3),'show')
%                 ylabel(ax_f(3),'Lac_P')
%                 title(strcat('Pyr_a_s_t=',pyr_ast(j)))
%                 end
%                 figure(f7)
%                 lgn = strcat('Lac_J=',num2str(v(3)/init_val(4)),'Pyr_A_0');
%                 hold(ax_f,'on')
%                 plot(t,flows(7,:),'DisplayName',lgn)
            end
            legnds = strcat('Pyr_A=',num2str(v(2)/init_val(2)),'Pyr_A_0');
            figure(f6)
            hold (ax_j,'on')
            plot(ax_j,lac_art,reshape(demand(i,j,:),length(lac_art),1),'LineWidth',1.5,'DisplayName',legnds)
            legend (ax_j,'show')
            xlabel('artery lactate (mM)')
            ylabel('pyramidal neurons dip (mM)')
            grid on
        end
%         figure
%         surf(pyr_ast,lac_art,reshape(demand(i,:,:),length(pyr_ast),length(lac_art)))
%         xlabel('pyr_ast')
%         ylabel('lac_art')
%         title('pyr_neuron')
%         zlabel('dip')
%         grid on
    end
    params.dyn.dip = demand;
end
%% Initial dip with different Vm_AE
if (multi == 4)% dynamics vs lacj vs V_m_AE
    c = 4;
    f5 = figure('Name','inputs');
    tau = params.dyn.tau;
    f6 = figure('Name','demand vs artery lac.','Position', [10 10 1000 600]);
    ax_j = subplot(1,1,1);
%     colors = [colors(1,:);colors(2,:);[1 1 1];[0 0.5 0.1]];
%     set(gca,'ColorOrder',colors)
    lac_art = linspace(params.Lac_j,params.Lac_j*10,c);
    VAE_vector = params.glb.MCT4.VAEs(1:3:20);
%     colors = colormap(winter(length(VAE_vector)));
%     colors = [zeros(1,length(VAE_vector));linspace(0,1,length(VAE_vector));linspace(1,0.5,length(VAE_vector))]';
    for i = 1:length(VAE_vector)
        params.glb.MCT4.fig_type = "typical";
        params.MCT.AE.Vm = VAE_vector(i);
        params = InitSS(params);
        params.pyr_P = params.glb.bestSS(8)/params.LDH.N1p;
        params.pyr_A = params.glb.bestSS(9)/params.LDH.ka;
        params.dyn.tau = [0.4 0.4 0.4 0.4 1 1 0.4 0.4];
        params.dyn.mode = 2;
        params.dyn.v = [params.pyr_P*0.98,params.pyr_A,params.CBF.F0*0.6,...
                        params.Lac_j];
        params.dyn.func = ["SHS","SHS","DE","SHS"];
        ss = [params.glb.bestSS(8),params.glb.bestSS(5),params.glb.bestSS(9),...
        params.glb.bestSS(6)];
        for k = 1:length(lac_art)
            v(1) = params.dyn.v(1);
            v(2) = params.pyr_A;
            v(3) = params.CBF.F0*1.5;
            v(4) = lac_art(k);
            [t,X] = ode15s(@(t,x)sys4D(t,x,params,0,[tau v]),params.ft,ss);
            flows = sys4D(t',X',params,1,[tau v],f5);
            if ((ss(1) - min(X(:,1)))>=0)
                demand(i,k) = (max(X(:,1)) - min(X(:,1)))/ss(1);
%                     demand_flows (i,j,k) = flows(1,1) - min() 
            else
                demand(i,k) = 0;
            end
%                 if (k==1)
            input =  tdf(params,[init_val(4),v(4),tau(end-1:end)],'DE2');
            lgn = strcat('Lac_J=',num2str(max(input)/params.Lac_j),'Lac_J_0');
        end
        legnds = strcat('V_A_E^m=',num2str(params.MCT.AE.Vm));
        figure(f6)
        hold (ax_j,'on')
        lacj_plot = plot(lac_art,demand(i,:),'LineWidth',1.5,'DisplayName',legnds);
        legend (ax_j,'show','Location','eastoutside')
        xlabel('artery lactate (mM)')
        ylabel('contribution of Lac_J on Lac_P')
%         h = get(gca, 'Children');
%         set(h(1), 'Color', colors(i,:));
%         set(lacj_plot,'Color',colors(i,:))
     end
end
box on
% params.dyn.dip = demand;

%% ------------------Visualization-----------------%
    function a1 = plotting(varargin)
        t1 = params.CBF.t1;
        tend = params.CBF.tend;
        dt = tend-t1;
        switch char(varargin)
            case 'flow'
                for r = 1:6
                   figure(f3) 
                   hold (axe(r),'on')
                   box(axe(r),'on')
                   if (r<=2)
                     a1(m,r) = plot(axe(r),t,flows(r,:));
%                     legend (axe(r),'show')
                   else
                       if (r == 3)
                           plot(axe(r),t,flows(r,:)-flows(r+1,:))
%                           legend (axe(r),'show')
                       end
                       if (r == 4)
                           plot(axe(r),t,flows(r+1,:))
%                           legend (axe(r),'show')
                       end
                       if (r>4)
                           hold (axe(r),'on')
                           plot(axe(r),t,flows(r+1,:))
                       end
                   end
                    if(m == c)
                        ylim(axe(r), 'padded')
                        g_x = t1:dt:tend; % user defined grid X [start:spaces:end]
                        g_y = get(axe(r),'ylim');
                        for l=1:length(g_x)
                            plot(axe(r),[g_x(l) g_x(l)],g_y,'k:') %y grid lines  
                            hold on
                            plot(axe(r),[t1 tend],[g_y(2) g_y(2)],'k:')
                        end
    %                     if(m==5)
    %                         plot(axe(r),[3 6],[g_y(2) g_y(2)],'r:')
    %                     end

%                         xticks(axe(r),([t1:dt:params.L]));
%                         xticklabels(axe(r),string([0:dt:params.L-dt]))
                        xticks(axe(r),[t1:dt:tend (tend+2*dt:dt*2:params.L)]);
                        xticklabels(axe(r),['0',string(dt),string(tend+2*dt-t1:dt*2:params.L-t1)]);
                   end
                end
                for jj=2:2
                    hold(axe(jj),'on')
                    legend(legnds) 
                end
            case 'LDH'
                figure(f4)
                    for r = 1:2
                        hold (ax4(r),'on')
                        plot(ax4(r),t,flows(r+5,:))
%                         legend (ax4(r),'show')
%                         grid(ax4(r),'on')
%                         if(m == c)
                            g_x = t1:dt:tend; % user defined grid X [start:spaces:end]
                            g_y = get(ax4(r),'ylim');
                            for l=1:length(g_x)
                                plot(ax4(r),[g_x(l) g_x(l)],g_y,'k:') %y grid lines  
                                hold on
                                plot(ax4(r),[t1 tend],[g_y(2) g_y(2)],'k:')
                            end
                            xticks(ax4(r),([t1:dt:params.L]));
                            xticklabels(ax4(r),string([0:dt:params.L-dt]))
%  
%                         end
                    end
                    for jj=2:2
                        hold(ax4(jj),'on')
                        legend(legnds) 
                    end
            case 'conc'
                for r = 1:4
                    figure(f1)
                    hold (ax(r),'on')
                    box(ax(r),'on')
                    a1(m,r) = plot(ax(r),t,X(:,r));
                    if(m == c)
                        ylim(ax(r), 'padded')
                        g_x = t1:dt:tend; % user defined grid X [start:spaces:end]
                        g_y = get(ax(r),'ylim');
                        for l=1:length(g_x)
                            plot(ax(r),[g_x(l) g_x(l)],g_y,'k:') %y grid lines  
                            hold on
                            plot(ax(r),[t1 tend],[g_y(2) g_y(2)],'k:')
                        end
%                         xticks(ax(r),([t1 tend]));
%                         xticklabels(ax(r),{'0','3'})
                        xticks(ax(r),[t1:dt:tend (tend+2*dt:dt*2:params.L)]);
                        xticklabels(ax(r),['0',string(dt),string(tend+2*dt-t1:dt*2:params.L-t1)]);
                    end
                end
                for kk=2:2
                    hold(ax(kk),'on')
                    legend(legnds) 
                end
        end
    end
    function decay_time = decay_calc(vec) 
        [val,ind]= max(vec);
        tpeak = t(ind,1);
        diff = val - 0.62*(val - 1);
        [~ ,closestI] = min(abs(vec(ind:end,1)-diff));
        decay_time = t(ind+closestI-1,1)-tpeak;
        params.sim.tau_c = decay_time;
        disp(decay_time)
    end
end