%% Initialization+SteadyState+loading data: optimization, visulization and selection %%
clear
close all
clc
params.glb.optSS = 0; %1:Optimization, 0: Load previous results
params.dyn.mode = 0;
params.glb.numofsim = 1500;
params.glb.method = "ICA";
params.glb.MCT4.fig_type = "MCT4";%"typical";%
% n = 10;
% VAE_vector = [5.6,7.6];

% for i=1:2
%     params.MCT.AE.Vm = VAE_vector(i)
    params = InitSS(params);
% end
%  sols(i,:,:,:) = params.glb.sols;
%%
close all
params.dyn.artery.mode = 0;
param_dist_VAE(params,1);


%% Dynamics__perturbed parameters
    close all
%     idx = 20;
%     params.glb.bestSS = params.glb.MCT4.vals{1,idx}.solica(1,:);127
    idx = 20;
    params.glb.bestSS = params.glb.MCT4.vals{1,idx}.solica(127,:);
    params.MCT.AE.Vm = params.glb.MCT4.VAE(idx);
    params.L = 24;
    params.dt = 0.0005;
    params.CBF.t1 = 3;
    params.CBF.tend = 6;% min: 3+1.6
    params.ft = 0:params.dt:params.L-params.dt;
    params.CBF.rep = 1;
    params.dyn.tau = [0.4 0.4 0.4 0.4 1 1 0.4 0.4];
    params.pyr_P = params.glb.bestSS(8)/params.LDH.N1p;
    params.pyr_A = params.glb.bestSS(9)/params.LDH.ka;
    params.dyn.mode = 2;
    params.dyn.v = [params.pyr_P*0.8,params.pyr_A,params.CBF.F0*0.6,...
                    params.Lac_j];%with pyr_A*8 we will have the dip clearer
    params.dyn.func = ["SHS","SHS","DE","SHS"];
    params = dynamics_sys4D(params,1);
