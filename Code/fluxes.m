function F = fluxes(params,z,varargin)
set(0,'defaultAxesFontSize',14,'defaultAxesFontName','Times',...
    'defaultAxesFontSmoothing','on',...
    'defaultLineLineWidth',1.3,'defaultAxesLineWidth',1.3,'defaultFigureColor','w')
MCT = params.MCT;
if (params.dyn.mode == 0)
    xp = z(8);
    xe = z(5);
    xa = z(9);
    xc = z(6);
    pyr_P = xp/params.LDH.N1p;
    pyr_A = xa/params.LDH.ka;
    bol = params.Lac_j;
    cbf = params.CBF.F0;
    Pa = z(7)*pyr_A./(0.084+pyr_A);
else 
    if (params.dyn.mode == 2 || params.dyn.mode == 3) 
        % Creating inputs
        cbf = tdf(params,params.dyn.cbf,params.dyn.func(3));
        bol = tdf(params,params.dyn.bol,params.dyn.func(4));
%         params.CBF.t1 = params.CBF.t1+0.2;
        pyr_P = tdf(params,params.dyn.pyr_p,params.dyn.func(1));
        pyr_A = tdf(params,params.dyn.pyr_a,params.dyn.func(2));
        % interpolation
        t = varargin{1,1};
        cbf = interp1(params.ft,cbf,t);
        bol = interp1(params.ft,bol,t);

        pyr_P = interp1(params.ft,pyr_P,t);
        pyr_A = interp1(params.ft,pyr_A,t);

        if (max(bol) == params.Lac_j)
            legnds = strcat('Lac_J=',num2str(ceil(min(bol)/params.Lac_j)),'Lac_J_0');
        else
            legnds = strcat('Lac_J=',num2str(ceil(max(bol)/params.Lac_j)),'Lac_J_0');
        end
        if (length(t)>1)
            f5 = varargin{1,3}{1,1};
            t1 = params.CBF.t1;
            tend = params.CBF.tend;
            dt = tend-t1;
            figure(f5)
            Names = ["Pyr_N (mM)";"Pyr_A (mM)";"CBF (min^-^1)";...
                    "Lac_J (mM)"];
            for a = 1:4
               ax(a) = subplot(1,4,a);
               ax(a).YLabel.FontWeight = 'bold';
               ax(a).YLabel.String = Names(a);
%                 ax(a).XLabel.String = 'time (min)';
            end
%             set(ax,'FontSize',16,'FontName','Calibri',...
%                 'FontSmoothing','on')% Fontweight: bold
            params.dyn.inputs = [pyr_P;pyr_A;cbf;bol];
            for i=1:4
               hold (ax(i),'on')
               plot(ax(i),t,params.dyn.inputs(i,:),'DisplayName',legnds,'LineWidth',1.3)
               box  (ax(i),'on')
               xticks(ax(i),[t1:dt:tend (tend+2*dt:dt*2:params.L)]);
               xticklabels(ax(i),['0',string(dt),string(tend+2*dt-t1:dt*2:params.L-t1)]);
               ylim(ax(i), 'padded')
%                if (i==2)
                    ylim(ax(i),[min(params.dyn.inputs(i,:))*0.9 ...
                        max(params.dyn.inputs(i,:))*1.1])
%                end
%                g_x = t1:dt:tend; % user defined grid X [start:spaces:end]
%                g_y = get(ax(i),'ylim');
%                for l=1:length(g_x)
%                    plot(ax(i),[g_x(l) g_x(l)],g_y,'k:') %y grid lines  
%                    hold on
%                    plot(ax(i),[t1 tend],[g_y(2) g_y(2)],'k:')
%                end
            end
        end
        % assigning variables
        x = varargin{1,2};
        if params.dyn.mode == 2
            xp = x(1,:);
            xe = x(2,:);
            xa = x(3,:);
            xc = x(4,:);
        end
        Pa = z(5)*pyr_A./(0.084+pyr_A);
    end
end
% LDH
Pp = z(1)*pyr_P./(pyr_P+z(2));
Cp = z(3)*xp./(z(4)+xp);
Ca = params.LDH.ast.Vn*xa./(xa+params.LDH.ast.Kn);
% MCT
V_ep = MCT.EP.Vm*(xe-xp)./(MCT.EP.Km+xp+xe);
V_ae = MCT.AE.Vm*(xa-xe)./(MCT.AE.Km+xe+xa);
V_ac = MCT.AB.Vm*xa./(MCT.AB.Km+xa);
V_ca = MCT.BA.Vm*xc./(MCT.BA.Km+xc);
V_cap = 2*cbf.*(bol - xc)/params.vol.VC;
% output
F = [V_ep;V_ae;V_ac;V_ca;V_cap;Pp;Cp;Pa;Ca];
end