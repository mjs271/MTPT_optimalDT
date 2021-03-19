%==========================================================================
% This script conducts an error analysis of the MTPT method over a range of 
% particle numbers (N) and time step lengths (dt) for a 1D problem with
% constant velocity.

% This script, will generate Figures 2 and 3 in:
%     "Optimal time step length for Lagrangian, interacting-particle simulations
%      of diffusive mixing," XXX Journal XXXX year.
%==========================================================================

%%

clear variables

%% set default figure settings and get the color/line/marker order for plotting

[figSettings, color, rgbcmy, mList, lineList, dashList] = setFigureDefaults();

%% simulation parameters

% number of dt values to loop over
nDt = 9;
% small and large (log) values for dt
smallDt = -5;
bigDt = 1;

% vector of particle numbers to loop over
nVec = [250 500 1000 2500 5000];
% % Add 1 particle to all for the Heaviside case
% nVec = nVec + 1;

% vector of dt's to loop over
dtVec = logspace(smallDt, bigDt, nDt);

% length of domain
L = 50;
% diffusion coefficient
D = 1e0;
% max simulation time
maxT = 10;

% beta from eq. (6)--controls bandwidth of kernel
beta = 1/2;

% optimal dt calculation
optDt = (L ./ nVec).^2 ./ ((1 / beta) * 2 * D);

% concatenate the optimal dt's with the logspaced dtVec
dtVec = sort([optDt dtVec]);

% length of the N and dt vectors
nNum = length(nVec);
nDt = length(dtVec);

% adjust dtVec so it evenly divides maxT
dtVec = maxT ./ round(maxT ./ dtVec);

% numbers of standard deviations to search for nearest neighbors with
% rangesearch
numSD = 6;

% variance for Gaussian IC
sigma = 1.5e-3 * sqrt(4 * D);

% cell array for holding errors for plotting
errVec = cell(nNum, nDt);

%% run the simulation

% particle number (N) loop
for idxN = 1 : nNum

%     print info to screen
    fprintf('N loop: %i of %i\n', idxN, nNum)

%     set N
    N = nVec(idxN);
%     space the particles evenly
    X = linspace(0, L, N)';

%     dt loop
    for idxDt = 1 : nDt

%         print info to screen
        fprintf('dt loop: %i of %i\n', idxDt, nDt)

%         set dt and calculate the number of steps to take
        dt = dtVec(idxDt);
        nSteps = maxT / dt;

%         Heaviside IC (assumes domain is [0, L])
        mass = zeros(N, 1);
        mass(X < 0.5 * L) = 1;
        analytic = 0.5 * erfc((X - L / 2) ./ sqrt(4 * D * maxT));
%         normalize both initial mass and analytic solution to sum to 1
        mass = mass / sum(mass);
        analytic = analytic / sum(analytic);

% %         gaussian IC (assumes domain is [0, L])
%         mass = (1 / sqrt(2 * pi * sigma^2)) * exp(-((0.5 * L - X).^2 / (2 * sigma^2)));
%         initial = mass;
%         analytic = (1 / sqrt(2 * pi * (sigma^2 + 2 * D * maxT))) *...
%                      exp(-((0.5 * L - X).^2 / (2 * (sigma^2 + 2 * D * maxT))));
% %         normalize both initial mass and analytic solution to sum to 1
%         mass = mass / sum(mass);
%         analytic = analytic / sum(analytic);

%         calculate the max interaction distance for rangesearch()
        dist = numSD * sqrt((1 / beta) * 2 * D * dt);

%         conduct the rangesearch to find nearby particles
        [idx, r] = rangesearch(X, X, dist, 'BucketSize', ceil(1e-2 * N));

%         determine how many particles are nearby and preallocate the vectors
%         to build the sparse weight matrix
        Nclose = sum(cellfun('length', idx));
        row = zeros(Nclose, 1);
        col = zeros(Nclose, 1);
        val = zeros(Nclose, 1);

%         calculate the entries of the weight matrix
        start = 1;
        for i = 1 : N
            finish = start - 1 + length(idx{i});
%             evaluate mass-transfer kernel
            row(start : finish) = idx{i};
            col(start : finish) = i;
            val(start : finish) = (1 / sqrt((1 / beta) * 4 * pi * D * dt)) .*...
                                  exp(-((r{i}).^2 ./ ((1 / beta) * 4 * D * dt)));
            start = finish + 1;
        end

%         when N is super large the sparse methods and clearing large vectors
%         saves memory
        clear idx r

%         create the sparse weight matrix
        Wmat = sparse(row, col, val);
        clear row col val

%         normalize Wmat such that it is symmetric, using the arithmetic mean of
%         colsum and rowsum--this method looks super goofy but is necessary to
%         keep things sparse and save memory
        colsum = sparse(sum(Wmat));
        rowsum = sparse(sum(Wmat, 2));
        Wmat = spfun(@(x) 1 ./ x, (2 .* Wmat) ./ rowsum) +...
               spfun(@(x) 1 ./ x, (2 .* Wmat) ./ colsum);
        clear rowsum colsum
        Wmat = spfun(@(x) 1 ./ x, Wmat);

%         build the transfer matrix
        Tmat = speye(N) - spdiags(sum(beta .* Wmat, 2), 0, N, N) + beta .* Wmat;

        clear Wmat

%         conduct the mass transfers
        for i = 1 : nSteps
            mass = Tmat * mass;
        end

%         save the error
        errVec{idxN, idxDt} = sqrt(mean((mass - analytic).^2));

    end
end

%% plotting

Mlist={'d', 'o', 's', 'x', '*'};
Nlegend = strcat('\boldmath $N = $', {' \bf'}, split(num2str(nVec)));

fig = 3;
figure(fig)
clf
hold on
box on
h = zeros(1, nNum);
dtC = zeros(1, nNum);
for i = 1 : nNum
    h(i) = plot(dtVec, [errVec{i, :}],  'color', color(i, :), 'marker', Mlist{i});
    dtC(i) = line([optDt(i) optDt(i)], [min([errVec{:}]) ...
         max([errVec{:}])], 'linestyle', '--',  'color', color(i, :));
end
fakeDash = plot(nan, nan, 'k--');
legend([h fakeDash], [Nlegend ; '\boldmath $\widehat{\Delta t} = \frac{(L / N)^2}{2 D \beta^{-1}}$'], 'location', 'southwest')
xlabel('\boldmath$\Delta t$')
ylabel('\textbf{RMSE}')
set(gca,'XDir','reverse','XScale','log','YScale','log');
axis tight
