clearvars
%parameters
seed = sum(fix(clock)); % random seed
rng(seed);
tmax = 10; % end time
ntype = 4;
population = [1, 0, 0, 0]; % initial population
alpha = [1, 2, 3, 0]; % birth rate
beta = [0, 0, 0, 0]; % death rate
u = 1E-3;
U = [0, 1, 1, 0; 0, 0, 0, 1;0, 0, 0, 1;0, 0, 0, 0]*u;
transition_rates = [alpha' beta' U];
events = zeros(ntype, ntype, ntype + 2);
events(:,:,1) = eye(4); % birth events
events(:,:,2) = -eye(4);
[r,c] = find(transition_rates(:,3:end) ~= 0); % find possible events
for i = 1:length(r)
    events(r(i),c(i),c(i)+2) = 1;
end
% record grid
record_time = 0:0.1:tmax;
traj = [record_time' zeros(length(record_time),ntype)];
t = 0; % starting time
record = 1;
jumps = 0;
% Gillespie Algorithm 
% for a discrete-state continous-time Markov Chain
while t < tmax
    if record_time(record) <= t
        traj(record,1) =t;
        traj(record,2:end) = population;
        record = record + 1;
    end
    weights = transition_rates .* population';
    [R,C] = Offspring(weights);
    [tau] = Lifespan(weights);
    t = t + tau;
    population = population + events(R,:,C);
    jumps = jumps + 1;
end
traj(record,1) = t;
traj(record,2:end) = population;
record = record + 1;

% Plot a traj
semilogy(traj(:,1),traj(:,2:end));
legend('type-0','type-1','type-2','type-3');
xlabel('time');
ylabel('population');
xlim([0,15]);

function [R,C] = Offspring(weights) % return the row and col of the event
    [r,c] = find(weights(:,:) ~= 0); % find possible one step events
    prob = weights(weights ~=0)';
    prob = prob/sum(prob);
    uniform_rv = rand;
    for i = 1:length(r)
        uniform_rv = uniform_rv - prob(i);
        if uniform_rv < 0
            break;
        end
    end
    R = r(i); C = c(i);
end

function [tau] = Lifespan(weights)
    m = 1/sum(sum(weights));
    tau = exprnd(m);
end
