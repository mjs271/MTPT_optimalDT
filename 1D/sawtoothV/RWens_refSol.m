clear variables

% define dt
dtRW = 1e-3;

% number of bins for generating concentrations
nBins = 51;

% number of random walk particles
nRW = 1e6;

% number of ensemble runs to conduct
nRWens = 5e2;

% define the travel time across a single velocity node
tNode = 2;

% diffusion coefficient
D = 1e-3;

% initial position of point source/interface
X0 = 2;

% number of velocity fields to  to loop over
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

    % total simulation time (depends on velocity field)
    maxTvec(i) = 2 / vMat(i, 1) + 2 * nVperiods * tNode;
end

% length of 1D domain (make large enough so particles don't reach boundary)
L = 2 * X0 + 6 * nVperiods * Lv;

% x-values for piecewise velocity grid function, defined below
xV = (0 : Lv : L)';

%% run the simulation

% arrays for holding quantities of interest
RWSolEns = zeros(nBins - 1, nV);
bins = zeros(nBins, nV);
LInit = zeros(1, nV);

for idxV = 1 : nV

%     piecewise velocity function
    vGrid = vMat(idxV, 1) * ones(size(xV));
    vGrid(2 : 2 : end) = vMat(idxV, 2);
    vGrid(1 : 4) = vMat(idxV, 1);
    vFunc = @(x) interp1(xV, vGrid, x, 'linear');

%===============================================================================
%                   RWPT SIMULATION (once per velocity field)
%===============================================================================

%     print some info to screen
    fprintf('===> Starting RWPT Simulation for Vel #%i\n', idxV)

%     ensemble loop
    for ens = 1 : nRWens

%         start the timer
        tic

%         adjust dtVec so it evenly divides maxT (which may be different for
%         each velocity field)
        maxT = maxTvec(idxV);
        dt = dtRW;
        dt = maxT ./ round(maxT ./ dt);

        nSteps = maxT / dt;

%         point-source IC (assumes domain is [0, L])
        xRW = X0 * ones(nRW, 1);

%         time stepping loop
        for idxTime = 1 : nSteps

%===============================================================================
%                   ADVECTION (RK3--KUTTA'S METHOD)
%===============================================================================

            k1 = vFunc(xRW);
            k1X = xRW + (dt / 2) * k1;
            k2 = vFunc(k1X);
            k2X = xRW + dt * (2 * k2 - k1);
            k3 = vFunc(k2X);
            xRW = xRW +  (dt / 6) * (k1 + 4 * k2 + k3);

%===============================================================================
%                   RANDOM WALK
%===============================================================================

            xRW = xRW + randn(nRW, 1) * sqrt(2 * D * dt);

        end

            if ens == 1
%                 get the length of the first RWPT plume, and create the
%                 corresponding bins
                LBin = max(xRW) - min(xRW);
                bins(:, idxV) = linspace(min(xRW), max(xRW), nBins);
            end

%             bin the particles to get concentrations
            [RWBin, ~, ~] = histcounts(xRW(:, 1), bins(:, idxV), 'Normalization', 'pdf');

            RWSolEns(:, idxV) = RWSolEns(:, idxV) + RWBin';

%             stop the timer and print info to screen
            timer = toc;
            fprintf('RWPT time (ens #%i) = %3.2f seconds \n', ens, timer)
    end

%     generate concentrations (assumes each particle has mass 1/N)
    RWSolEns(:, idxV) = RWSolEns(:, idxV) / nRWens;

%===============================================================================

%     initial length of particle injection (calibrate to final width of RWPT plume)
    LInit(idxV) = LBin;

end

%% save the quantities of interest

RWref = struct('bins', bins, 'RWSols', RWSolEns, 'initLength', LInit);
save('RWrefStruct.mat', 'RWref')
