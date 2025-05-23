function [Best, fBest, history] = CEO(func, Np, Dim, Varmin, Varmax, N, MaxFES)
rand('state', sum(100*clock));

if mod(Np,2)~=0
    error('Np must be even！')
end

% Initialize the main population
Population = repmat(Varmin,Np, 1) + rand(Np, Dim) .* (repmat(Varmax - Varmin, Np, 1));
fit = fitness(func,Population); 
[fBest,index_best] = min(fit);
Best = Population(index_best,:);
history(1) = fBest;

% chaotic initial search domain
low_chacos = [-0.5 -0.25]; 
up_chacos = [0.5 0.25];
FEvals = Np;
t = 1; 
while FEvals < MaxFES
    oldfBest = fBest;
    rand_num = randperm(Np);
    ub = max(Population); % Upper limit of population per iteration
    lb = min(Population); % Lower limit of population per iteration
    
    for i = 1:2:Np
        % Randomly select two individuals
        index = rand_num(i:i+1);
        xy = Population(index,:); 
        
        % Perform interval mapping on xt and yt by executing Eq. (4).
        xy_dot = (xy -  lb)./((ub - lb)+eps).*(repmat(up_chacos',1,Dim) - repmat(low_chacos',1,Dim)) + repmat(low_chacos',1,Dim);
        
        % N chaotic individuals are obtained by executing Eq. (2)
        chaos_total = EDM(xy_dot(1,:),xy_dot(2,:),N);
        
        % Evaluate each chaotic sequence
        for k = 1 : 2 
            xy_chaos = chaos_total((k-1)*N+1 : k*N,:);
            % Executing Eq. (5) yields the actual position of the corresponding optimization problem.
            xy_chaos_dot = (( xy_chaos -  repmat(low_chacos(:,k),1,Dim) )./(repmat(up_chacos(:,k),1,Dim) - repmat(low_chacos(:,k),1,Dim) )).*(ub - lb) + lb; 
            xy_chaos_dot = boundConstraint (xy_chaos_dot, Varmin, Varmax);
            
            if rand < 0.5 % 
                xy_hat = xy(k,:) + rand(N,1).*( xy_chaos_dot - xy(k,:) ); % Mutation Eq. (7)
            else
                xy_hat = Best + rand(N,1).*( xy_chaos_dot - xy(k,:) ); % Mutation Eq. (8)
            end
            
            CR = rand;
            xy_trial = Binomial_crossover(xy(k,:), xy_hat, CR); % Crossover Eq. (9)
            xy_trial = boundConstraint (xy_trial, Varmin, Varmax);
            fit_xy_trial = fitness(func,xy_trial);
            [fBest_xy_trial,index_best] = min(fit_xy_trial);
            xy_trial_star = xy_trial(index_best,:);
            if  fBest_xy_trial < fit(index(k)) % Selection Eq. (10)
                Population(index(k),:) = xy_trial_star;
                fit(index(k)) = fBest_xy_trial;
            end   
        end
        
    end
    
    [fBest,index_best] = min(fit);
    Best = Population(index_best,:);

    if handleStagnation(oldfBest, fBest)
        break; 
    end
    
    t = t + 1;
    if mod(t,10)==0
        fprintf('iter=%d  ObjVal=%20.16g\n',t,fBest);
    end
    history(t) = fBest;
    
    FEvals = FEvals + N*Np;
end

end

function xy = EDM(x0, y0, iter)
% exponential discrete memristor (E-DM) map
k = 2.66;
x(1,:) = x0;
y(1,:) = y0;
for j = 1:iter
    x(j+1,:) = k * (exp(-cos(pi*y(j,:))) - 1) .* x(j,:);
    y(j+1,:) = y(j,:) + x(j,:);
end
x = x(2:end,:); y = y(2:end,:);
xy = [x;y]; % Merging particles created by chaos
end

function u = Binomial_crossover(p, v, CR)
% Binomial crossover
[N,dim] = size(v);
t = rand(N, dim) < CR;
random_cols = randi(dim, [N, 1]);
linear_indices = sub2ind([N, dim], (1:N)', random_cols);
t(linear_indices) = 1;
u = t .* v + (1 - t) .* p;
end

function fit = fitness(func,x)
% calculate fitness
N = size(x,1);
fit = zeros(N,1);
for i = 1:N
    fit(i) = feval(func,x(i,:));
end
end

function v = boundConstraint (v, Varmin, Varmax)
% Handle the elements of the state vector which violate the boundary
low_state = repmat(Varmin,size(v,1),1);
up_state = repmat(Varmax,size(v,1),1);
vioLow = v < low_state;
v(vioLow) = min(up_state(vioLow), 2 .* low_state(vioLow) - v(vioLow));
vioUpper = v > up_state;
v(vioUpper) = max(low_state(vioUpper), 2 .* up_state(vioUpper) - v(vioUpper));
end

function stop = handleStagnation(old, new)
persistent counter
if isempty(counter), counter = 0; end
if abs(old - new) < 1e-8
    counter = counter + 1;
    stop = (counter > 50);
else
    counter = 0;
    stop = false;
end
end

