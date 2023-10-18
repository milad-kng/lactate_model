%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA118
% Project Title: Implementation of Imperialist Competitive Algorithm (ICA)
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

function [BestSol, BestCost] = ica(fun,nv,lb,ub)
%% Problem Definition

CostFunction=fun;        % Cost Function

nVar=nv;             % Number of Decision Variables

VarSize=[1 nVar];   % Decision Variables Matrix Size

VarMin=lb;         % Lower Bound of Variables
VarMax= ub;       % Upper Bound of Variables


%% ICA Parameters

MaxIt=2000;         % Maximum Number of Iterations

nPop=50;            % Population Size
nEmp=15;            % Number of Empires/Imperialists

alpha=1;            % Selection Pressure

beta=1.5;           % Assimilation Coefficient

pRevolution=0.05;   % Revolution Probability
mu=0.1;             % Revolution Rate

zeta=0.2;           % Colonies Mean Cost Coefficient


%% Globalization of Parameters and Settings

global ProblemSettings;
ProblemSettings.CostFunction=CostFunction;
ProblemSettings.nVar=nVar;
ProblemSettings.VarSize=VarSize;
ProblemSettings.VarMin=VarMin;
ProblemSettings.VarMax=VarMax;

global ICASettings;
ICASettings.MaxIt=MaxIt;
ICASettings.nPop=nPop;
ICASettings.nEmp=nEmp;
ICASettings.alpha=alpha;
ICASettings.beta=beta;
ICASettings.pRevolution=pRevolution;
ICASettings.mu=mu;
ICASettings.zeta=zeta;


%% Initialization

% Initialize Empires
emp=CreateInitialEmpires();
% load('saved/emp_march25.mat','emp')
% Array to Hold Best Cost Values
BestCost=zeros(MaxIt,1);


%% ICA Main Loop

for it=1:MaxIt
%     tic

    % Assimilation
    emp=AssimilateColonies(emp);
    
    % Revolution
    emp=DoRevolution(emp);
    
    % Intra-Empire Competition
    emp=IntraEmpireCompetition(emp);
    
    % Update Total Cost of Empires
    emp=UpdateTotalCost(emp);
    
    % Inter-Empire Competition
    emp=InterEmpireCompetition(emp);
    
    % Update Best Solution Ever Found
    imp=[emp.Imp];
    [~, BestImpIndex]=min([imp.Cost]);
    BestSol=imp(BestImpIndex);
%     save('saved/new/Final_Lac_P_art_pyr_P.mat','BestSol')
    % Update Best Cost
    BestCost(it)=BestSol.Cost;
%     save('saved/emp_march25.mat','emp')
%     save('saved/new/emp_LacPLacA.mat','emp') 
    % Show Iteration Information
%     disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
%     toc
end


%% Results

% figure;
% plot(sqrt(BestCost),'LineWidth',2);
% semilogy(sqrt(BestCost),'LineWidth',2);
% xlabel('Iteration');
% ylabel('Best Cost');
% grid on;
end