function params = glblopt(params,bounds,numofsim)

%% Create problem structure

lb = bounds(1,:);
ub = bounds(2,:);
nvars = length(lb);
numofmethods = 2;
%%
fval = zeros(numofsim,numofmethods);
xsol = zeros(numofsim,numofmethods,nvars);
dt = zeros(numofsim,numofmethods);
fval4 = zeros(numofsim,4);
multi = 0;
objfun = @(z)syslac(z,params,multi);% multi=0: sum of squares
multi = 1;
objfunmul = @(z)syslac(z,params,multi);
A = [1 0 0 0 0 0 -1 0 0];
b = 0;
for k = 1:numofsim
        %% Initialization
        z0 = lb + (ub-lb).*rand(1,nvars);
%         problem = createOptimProblem('fmincon','x0',z0,'objective',objfun,'lb',lb,'ub',ub,'nonlcon',@nonlc);
        sol = 1;
        %% Optimization
        %ICA
        addpath('ICA')
        t0 = tic();
        [bestsol, ~] = ica(objfun,length(z0),lb,ub);
        dt(k,sol) = toc(t0);
        fval(k,sol) = bestsol.Cost;
        disp(fval(k,sol))
        xsol(k,sol,:) = bestsol.Position;
        fprintf('done iteration number %i\n for ICA', k)
        delta = toc(t0);
        disp(delta)
        % fgoalattain
%         sol = sol+1;
%         t0 = tic();
%         goal = [0,0,0,0];
%         weight = [1,1,1,1];
%         options = optimoptions('fgoalattain','GoalsExactAchieve',4);
%         [xsol(k,sol,:), fval4(k,:),~] = fgoalattain(objfunmul,z0,goal,weight,A,b,[],[],lb,ub,@nonlc,options);
%         fprintf('done iteration number %i\n for FG', k)
        delta = toc(t0);
        disp(delta)
end
%     gsol = reshape(xsol(:,1,:),numofsim,nvars);
%     msol = reshape(xsol(:,2,:),numofsim,nvars);
%     isol = reshape(xsol(:,3,:),numofsim,nvars);
%     params.glb.timing = dt;
    params.glb.sols = xsol;
%     params.glb.solmean = [mean(gsol);mean(msol);mean(isol)];
%     params.glb.solstd = [std(gsol);std(msol);std(isol)];
    params.glb.fvalmean = mean(fval,1);
    params.glb.fvalstd = std(fval,1);
%     params.glb.solms = solutionms;
    params.glb.fval = fval;
%     params.glb.fval4 = fval4;
%     params = SelectSS(params);

end
