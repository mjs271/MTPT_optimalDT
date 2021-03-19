%==========================================================================
% This script conducts an error analysis of the MTPT method for a 1D problem
% with a "sawtooth" velocity field and considers 3 velocity pairs.

% Included in this directory are a script that calculates the v2 velocity
% and two .mat files containing reference solution data and experimentally
% derived optimal dt data.

% This script, will generate Figure 4 in:
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
smallDt = -3;
bigDt = 0;

% number of bins, so as to compare MT and RW results
nBins = 51;

% number of mass transfer and random walk particles
% (odd numbers to capture center of Gaussian)
nMT = 1e3 + 1;
nRW = 1e4;

% numbers of standard deviations to search for nearest neighbors with
% rangesearch
numSD = 6;

% beta from eq. (6)--controls bandwidth of kernel
beta = 1;

% vector of dt's to loop over
dtVec = logspace(smallDt, bigDt, nDt);

% define the travel time across a single velocity node
tNode = 2;

% diffusion coefficient
D = 1e-3;

% initial position of point source/interface
X0 = 2;

% number of velocity fields to loop over
nV = 3;
% number of velocity periods to cover in a simulation
nVperiods = 2;
% length of each velocity segment (note: not length of full velocity cycle--that would be 2 * Lv)
Lv = 1;
% define the v1's and calculate the v2's such that the travel time is equal to tNode
vMat = zeros(3, 2);
% these are the v1's
vMat(:, 1) = [0.25 0.75 1];

% total simulation time for a given velocity field
maxTvec = zeros(nV, 1);

for i = 1 : nV
%     calculate the v2's and organize so the ones in the first column are
%     always the smaller
    vMat(i, 2) = calcV2(vMat(i, 1), Lv, tNode);
    vMin = min(vMat(i, :));
    vMax = max(vMat(i, :));
    vMat(i, 1) = vMin;
    vMat(i, 2) = vMax;

    % calculate total simulation time (depends on velocity field)
    maxTvec(i) = 2 / vMat(i, 1) + 2 * nVperiods * tNode;
end

% length of 1D domain (make large enough so particles don't reach boundary)
L = 2 * X0 + 6 * nVperiods * Lv;

% cell array for holding errors for plotting
errVecMT = cell(nDt, nV);

% x-values for piecewise velocity grid function, defined below
xV = (0 : Lv : L)';

% load the RW results
load('RWrefStruct_N500.mat')
bins = RWref.bins;
RWSolEns = RWref.RWSols;
LInit = RWref.initLength;

% calculate the bins for plotting
plotBins = bins(2 : end, :) - (bins(2, :) - bins(1, :)) / 2;

% estimate optimal dt calculation to include in the simulation
% optDt = (LInit ./ nMT).^2 ./ ((1 / beta) * 2 * D);

% load experimentally-calculated optimal dt's (calculated below, but we
% want to include these values in the simulation)
load('expDs.mat')
optDtExp = expDs.beta1;

% concatenate the defined dtVec with optDt and redefine nDt
% dtVec = sort([optDtExp optDt dtVec]);
dtVec = sort([optDtExp dtVec]);
nDt = length(dtVec);

% track average and max inter-particle spacing for each simulation
% each element is nSteps x 2, and first column is average spacing, second
% is max spacing
dsVec = cell(nDt, nV);

%% run the simulation

% counter for printing simulation progress to screen
runCounter = 0;

% dtVec for a given velocity field
dtVecVel = zeros(length(dtVec), nV);

for idxV = 1 : nV

%     piecewise velocity function
    vGrid = vMat(idxV, 1) * ones(size(xV));
    vGrid(2 : 2 : end) = vMat(idxV, 2);
    vGrid(1 : 4) = vMat(idxV, 1);
    vFunc = @(x) interp1(xV, vGrid, x, 'linear');

%     adjust dtVec so it evenly divides maxT (which may be different for
%     each velocity field)
    maxT = maxTvec(idxV);
    dtVecVel(:, idxV) = maxT ./ round(maxT ./ dtVec);


%===============================================================================
%                   MTPT ENSEMBLE (over dt)
%===============================================================================

%     dt loop
    for idxDt = 1 : nDt

%         start timer for a given simulation and print some info to screen
        tic
        runCounter = runCounter + 1;
        fprintf('===> Starting MTPT Simulation for: Vel #%i, dt #%i (%i of %i)\n',...
                idxV, idxDt, runCounter, nV * nDt)

%         set the current dt
        dt = dtVecVel(idxDt, idxV);

%         caluclated the number of steps (guarantee that it's an integer)
        nSteps = round(maxT / dt);

%         each element is nSteps x 2, and first column is average spacing,
%         second is max spacing
        dsVec{idxDt, idxV} = zeros(nSteps, 2);

%         initial particle positions
        xMT = linspace(X0 - LInit(idxV) / 2, X0 + LInit(idxV) / 2, nMT)';

%         gaussian IC (assumes domain is [0, L])
        mass = zeros(nMT, 1);
        mass(ceil(nMT / 2)) = 1;
        initial = mass;

        curTime = 0;

%         time stepping loop
        for idxTime = 1 : nSteps

%===============================================================================
%                   ADVECTION (RK3--KUTTA'S METHOD)
%===============================================================================

            k1 = vFunc(xMT);
            k1X = xMT + (dt / 2) * k1;
            k2 = vFunc(k1X);
            k2X = xMT + dt * (2 * k2 - k1);
            k3 = vFunc(k2X);
            xMT = xMT +  (dt / 6) * (k1 + 4 * k2 + k3);

%===============================================================================
%                   MASS TRANSFER
%===============================================================================
%             calculate the max interaction distance for rangesearch()
            dist = numSD * sqrt((1 / beta) * 2 * D * dt);

%             conduct the rangesearch to find nearby particles
            [idx, r] = rangesearch(xMT, xMT, dist, 'BucketSize', ceil(1e-2 * nMT));

%             first column is average spacing, second is max spacing
            dsVec{idxDt, idxV}(idxTime, 1) = (max(xMT) - min(xMT)) / nMT;
            dsVec{idxDt, idxV}(idxTime, 2) = max(cell2mat(cellfun(@(x) min(x(x > 0)), r, 'uniformOutput', false)));

%             determine how many particles are nearby and preallocate the vectors
%             to build the sparse weight matrix
            Nclose = sum(cellfun('length', idx));
            row = zeros(Nclose, 1);
            col = zeros(Nclose, 1);
            val = zeros(Nclose, 1);

%             calculate the entries of the weight matrix
            start = 1;
            for i = 1 : nMT
                finish = start - 1 + length(idx{i});
%                 evaluate mass-transfer kernel
                row(start : finish) = idx{i};
                col(start : finish) = i;
                val(start : finish) = (1 / sqrt((1 / beta) * 4 * pi * D * dt)) .*...
                                      exp(-((r{i}).^2 ./ ((1 / beta) * 4 * D * dt)));
                start = finish + 1;
            end

%             when N is super large the sparse methods and clearing large vectors
%             when done saves memory
            clear idx r

%             create the sparse weight matrix
            Wmat = sparse(row, col, val);
            clear row col val

%             normalize Wmat such that it is symmetric, using the arithmetic mean of
%             colsum and rowsum--this looks super goofy but is necessary to
%             keep things sparse and save memory
            colsum = sparse(sum(Wmat));
            rowsum = sparse(sum(Wmat, 2));
            Wmat = spfun(@(x) 1 ./ x, (2 .* Wmat) ./ rowsum) +...
                   spfun(@(x) 1 ./ x, (2 .* Wmat) ./ colsum);
            clear rowsum colsum
            Wmat = spfun(@(x) 1 ./ x, Wmat);

%             build the transfer matrix
            Tmat = speye(nMT) - spdiags(sum(beta .* Wmat, 2), 0, nMT, nMT) + beta .* Wmat;

            clear Wmat

%             conduct the mass transfers
            mass = Tmat * mass;
        end

%         bin the masses from the MT simulation
        [~, ~, binMT] = histcounts(xMT, bins(:, idxV));
        MTSol = zeros(nBins - 1, 1);
        for i = 1 : nMT
            if binMT(i) == 0
                continue
            end
           MTSol(binMT(i)) = MTSol(binMT(i)) + mass(i);
        end

%         calculate concentrations for comparison
        MTSol = MTSol ./ (bins(2, idxV) - bins(1, idxV));

%         stop the timer and print run time
        timer = toc;
        fprintf('MTPT time = %3.2f seconds \n', timer)

%         normalize the solutions to sum to 1
        MTSol = MTSol / sum(MTSol);
        RWSolEns(:, idxV) = RWSolEns(:, idxV) / sum(RWSolEns(:, idxV));

%         calculate and save the error
        errVecMT{idxDt, idxV} = sqrt(mean((MTSol - RWSolEns(:, idxV)).^2));

    end
end

%% calculate optimal dt, based on the ds's calculated above

% optimal dt calculation based on particle spacings during run
maxdsVec = cellfun(@max, dsVec, 'UniformOutput', false);
maxdsVec = cell2mat(maxdsVec);
maxdsVec = max(maxdsVec);
optDtAvgDs = (maxdsVec(:, 1 : 2 : 5)).^2 ./ ((1 / beta) * 2 * D);

%% plotting

Mlist={'d', 'o', 's', 'x', '*'};
vLegend = {'\boldmath $V_1$', '\boldmath $V_2$', '\boldmath $V_3$'};
dtLegend = strcat('dt = ', {' '}, split(num2str(dtVec)));

fig = 3;
figure(fig)
clf
hold on
box on
h = zeros(1, nV);
dtC = zeros(1, nV);
dtMRun = zeros(1, nV);
for i = 1 : nV
    h(i) = plot(dtVec, [errVecMT{:, i}],  'color', color(i, :), 'marker', Mlist{i});
    dtMRun(i) = line([optDtAvgDs(i) optDtAvgDs(i)], [min([errVecMT{:}]) ...
                     max([errVecMT{:}])], 'linestyle', '--',  'color', color(i, :));
end
fakeDash = plot(nan, nan, 'k--');
legend([h fakeDash], [vLegend '\boldmath $\widehat{\Delta t} = \frac{\overline{\Delta s}^2_{\mathrm{max}}}{2 D \beta^{-1}}$'],...
       'location', 'northwest')
legend([h fakeDash], [vLegend '\boldmath $\widehat{\Delta t} = \frac{\overline{\Delta s}^2_{\mathrm{max}}}{2 D \beta^{-1}}$'],...
       'position', [0.18 0.8 0.1 0.1])
xlabel('\boldmath$\Delta t$')
ylabel('\textbf{RMSE}')
set(gca,'XDir', 'reverse', 'XScale', 'log', 'YScale', 'log');
axis tight
