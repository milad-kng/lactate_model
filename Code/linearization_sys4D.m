%% Linear Stability
% Jacobian
function [eigvec,eigval,W] = linearization_sys4D(params,X)

x = sym('x',[1 4]);
F = sys4D_csp (x,params,0,0);
%F = systest(x);
J = jacobian(F,x);
x1 = X(:,1);%1.5*10^-4;%ss(1);
x2 = X(:,2);%0.1*10^-6;%ss(2);
x3 = X(:,3);%300;%ss(3);
x4 = X(:,4);%ss(4);
jacsym = subs(J);
for i=1:3
    for j=1:3
        jacss(i,j) = sym2poly(jacsym(i,j));
    end
end
[eigvec,eigval,W] = eig(jacss);
%% Dynamic
% F = @(t,x) sys4D_nonD (x,params);
% [t,X]=NaiveRK4Sto(F,params.L,params.dt,0,ss);
% figure
% plot(t,X)