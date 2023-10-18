function x = tdf(params,varargin)

%HS: Heavyside, SHS: slope heavyside, DE: Double exponential
td = 0.2;
if (params.CBF.rep>1)
    L = params.CBF.tend;
else
    L = params.L;
end
dt = params.dt;
tdelay = 0.1;
% if (nargin == 3)
%     t1 = params.CBF.t1 + tdelay;
%     tend = params.CBF.tend + tdelay;
% else
    t1 = params.CBF.t1;
    tend = params.CBF.tend;
% end
t2 = t1 + td;
features = varargin{1,1};
func = varargin{1,2};
v0 = features(1);
v1 = features(2);

switch func
case 'comb_demand'
    x = v0*ones(L/dt,1);
    for i = 1:L/dt
        if((i <= t2/dt)&&(i > t1/dt))
            x(i) =  v0 + (v1-v0)*(i*dt-t1)/(t2-t1);
        else
            if (i>t2/dt)&&(i<tend/dt)
                x(i) = v1;
        else
            if (i>=tend/dt)&&(i<=(tend+td)/dt)
                x(i) = v1 + (-v1+v0)*(i*dt-tend)/(t2-t1);
            end
            end
        end
    end
    
case 'HS'
    x = v0*ones(L/dt,1);
    for i = 1:L/dt
       if ((i<tend/dt)&&(i>t1/dt)) 
           x(i,1) = v1;
       end
    end   
    
case  'SHS'
    x = v0*ones(ceil(L/dt),1);
    for i = 1:L/dt
        if((i <= t2/dt)&&(i > t1/dt))
            x(i) =  v0 + (v1-v0)*(i*dt-t1)/(t2-t1);
        else
            if (i>t2/dt)&&(i<tend/dt)
                x(i) = v1;
        else
            if (i>=tend/dt)&&(i<=(tend+td)/dt)
                x(i) = v1 + (-v1+v0)*(i*dt-tend)/(t2-t1);
            end
            end
        end
    end
    
case 'DE'
    x = v0*ones(L/dt,1);
%     tau(1) = features(3);
%     tau(2) = features(4);
%     tau(3) = features(5);
%     if(nargin == 7)
        for i = 1:L/dt
            if((i >= t1/dt)&&(i<tend/dt))
                x(i) =  v0 + v1*(exp((-i*dt+t1)/2) - exp((-i*dt+t1)/0.1));%2, 0.1
            else
                if (i>=tend/dt)
                    x(i) = v0*(1 + (x(ceil(tend/dt))-v0)*exp((tend-i*dt)/20));%20
                end
            end
        end
%     else
%         error("not enough input arguments")
%     end
    case 'DEN'
    x = v0*ones(L/dt,1);
%     tau(1) = features(3);
%     tau(2) = features(4);
%     tau(3) = features(5);
%     if(nargin == 7)
        for i = 1:L/dt
            if((i >= t1/dt)&&(i<tend/dt))
                x(i) =  v0*(1 - (v1/v0)*(exp((-i*dt+t1)/5) - exp((-i*dt+t1)/0.2)));
            else
                if (i>=tend/dt)
                    x(i) = v0*(1 + (x(ceil(tend/dt))-v0)*exp((tend-i*dt)/10));
                end
            end
        end
%     else
%         error("not enough input arguments")
%     end
    
case 'DE2'

    x = v0*ones(L/dt,1);
    tau(1) = features(3);
    tau(2) = features(4);
%     if(nargin == 6)
        for i = 1:L/dt
            if((i >= t1/dt)&&(i<=tend/dt))
                x(i) =  v0 +v1*(1-(exp((-(i*dt-t1))/tau(1))));%v1=9, 0.0059
            else
                if (i>tend/dt)
                    x(i) = v0 + (x(tend/dt)-v0)*exp((tend-i*dt)/tau(2));%1.0666
                end
            end
        end
%     else
%         error("not enough input arguments")
%     end

case 'DE2_exp'
    tau = params.exp.tau;
    x = v0*ones(L/dt,1);
    t2 = 8;
    tend = 18.5;
    for i = 1:L/dt
        
        if((i >= t2/dt)&&(i<=tend/dt))
            x(i) =  v0 +v1*(1-(exp((-(i*dt-t2))/tau(2))));
        else
            if (i>tend/dt)
                x(i) = v0 + (x(tend/dt)-v0)*exp((tend-i*dt)/tau(3));
            end
        end
    end
    
case 'DEN2'
    x = v0*ones(L/dt,1);
    tau(1) = features(3);
    tau(2) = features(4);
%     if(nargin == 6)
        for i = 1:L/dt
            if((i >= t1/dt)&&(i<=tend/dt))
                x(i) =  v0 -v1*(1-(exp((-(i*dt-t1))/tau(1))));
            else
                if (i>tend/dt)
                    x(i) = v0 + (x(i-1)-v0)*exp((tend-i*dt)/tau(2));
                end
            end
        end
%     else
%         error("not enough input arguments")
%     end
case 'DIP'
    x = v0*ones(L/dt,1);
    for i = 1:L/dt
        if (i>(t2-1)/dt)&&(i <= t2/dt)
            x(i) = v0*exp(((t2-1)-i*dt)/5);
            
        else
            if((i > t2/dt)&&(i<tend/dt))
%                 x(i) =  x(t0/dt) +v1*(1-(exp((-(i*dt-t0))/2)));
                x(i) =  v0 +v1*(1-(exp((-(i*dt-t2))/3)));
        else
            if (i>=tend/dt)
                x(i) = v0 + (x(tend/dt)-v0)*exp((tend-i*dt)/5);
            end
            end
        end
    end
end
end
    
