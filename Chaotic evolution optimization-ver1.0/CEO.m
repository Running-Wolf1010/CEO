function [Best, fBest, history] = CEO(func, Np, Dim, Varmin, Varmax, N)
rand('state', sum(100*clock));

if mod(Np,2)~=0
    error('Np must be set to an even number greater than 2！')
end

% Search Range
if length(Varmin)== 1
    lu = repmat([Varmin; Varmax], 1, Dim);
else
    lu = [Varmin; Varmax];
end

% Initialize the main population
[Population,fit,fBest,Best] = Initialization(func,lu, Np, Dim);
history(1) = fBest;

% chaotic initial search domain
low_chacos = [-0.5 -0.25]; 
up_chacos = [0.5 0.25];

FEvals = Np;
t = 1; % Initialization iteration number
counter = 0;
while 1
    oldfBest = fBest;
    rand_num = randperm(Np);
    ub = max(Population); % Upper limit of population per iteration
    lb = min(Population); % Lower limit of population per iteration
    
    for i = 1:2:Np
        index = rand_num(i:i+1);
        xy = Population(index,:); % Randomly select two individuals
        
        % Perform interval mapping on xt and yt by executing Eq. (4).
        xy_dot = ((xy -  lb)./(ub - lb)).*(repmat(up_chacos',1,Dim) - repmat(low_chacos',1,Dim)) + repmat(low_chacos',1,Dim);
        
        % N chaotic individuals are obtained by executing Eq. (2)
        [x_chaos,y_chaos] = EDM(xy_dot(1,:),xy_dot(2,:),N);
        
        chaos_total = [x_chaos;y_chaos]; % Merging particles created by chaos
        
        for k = 1 : 2 % Evaluate each chaotic sequence
            xy_chaos = chaos_total((k-1)*N+1 : k*N,:);
            % Executing Eq. (5) yields the actual position of the corresponding optimization problem.
            xy_chaos_dot = (( xy_chaos -  repmat(low_chacos(:,k),1,Dim) )./(repmat(up_chacos(:,k),1,Dim) - repmat(low_chacos(:,k),1,Dim) )).*(ub - lb) + lb; 
            
            if rand < 0.5 % 
                xy_hat = xy(k,:) + rand(N,1).*( xy_chaos_dot - xy(k,:) ); % Mutation Eq. (7)
            else
                xy_hat = Best + rand(N,1).*( xy_chaos_dot - xy(k,:) ); % Mutation Eq. (8)
            end
            
            CR = rand ;
            xy_trial = Binomial_crossover(xy(k,:), xy_hat, CR); % Crossover Eq. (9)
            xy_trial = boundConstraint (xy_trial, lu);
            fit_xy_trial = fitness(func,xy_trial);
            [fBest_xy_trial,index_best] = min(fit_xy_trial);
            xy_trial_star = xy_trial(index_best,:);
            if  fBest_xy_trial < fit(index(k)) % Selection Eq. (10)
                Population(index(k),:) = xy_trial_star;
                fit(index(k)) = fBest_xy_trial;
            end
            FEvals = FEvals + N;
        end
        
    end
    
    [fBest,index_best] = min(fit);
    Best = Population(index_best,:);
    
    t = t + 1;
    
    %  termination conditions
    if norm(oldfBest-fBest) < 1e-8 % can be changed
        counter = counter + 1;
        if counter > 50 % can be changed
            disp('满足停止条件');
            break;
        end
    else
        counter = 0;
    end
    
    if mod(t,100)==0
        fprintf('iter=%d  ObjVal=%g\n',t,fBest);
    end
    history(t) = fBest;
    
end

end

function [X,Y] = EDM(x0,y0,itermax)
% exponential discrete memristor (E-DM) map
k = 2.66;
x = x0; % trajectory domain [-1,1]
y = y0;% trajectory domain [-0.5,0.5]
xo = x;
yo = y;
% System iteration
for j = 1:itermax
    xn = k*(exp(-cos(pi.*yo))-1).*xo;
    yn = yo + xo;
    % Stored iteration value
    x =[x; xn];
    y =[y; yn];
    % Update the initial value of each iteration
    xo = xn;
    yo = yn;
end
X =x(2:end,:);
Y =y(2:end,:);
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
SE = size(x,1);
fit = zeros(SE,1);
for i = 1:SE
    fit(i) = feval(func,x(i,:));
end
end

function [population,fit,fBest,Best] = Initialization(func,lu, Np, Dim)
population = repmat(lu(1, :),Np, 1) + rand(Np, Dim) .* (repmat(lu(2, :) - lu(1, :), Np, 1)); % 初始化种群
fit = fitness(func,population); % Evaluate the objective function values
[fBest,index_best] = min(fit);
Best = population(index_best,:);
end

function v = boundConstraint (v, lu)
% Handle the elements of the state vector which violate the boundary
low_state = repmat(lu(1, :),size(v,1),1);
up_state = repmat(lu(2, :),size(v,1),1);
vioLow = v < low_state;
v(vioLow) = min(up_state(vioLow), 2 .* low_state(vioLow) - v(vioLow));
vioUpper = v > up_state;
v(vioUpper) = max(low_state(vioUpper), 2 .* up_state(vioUpper) - v(vioUpper));
end



