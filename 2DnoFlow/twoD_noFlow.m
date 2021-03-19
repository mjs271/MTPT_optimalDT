%==========================================================================
% This script conducts an error analysis of the MTPT method over a range of 
% particle numbers (N) and time step lengths (dt) for a 2D problem with
% constant velocity.

% This script, will generate Figures 5 and 6 in:
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
nVec = [51 75 101 125].^2;
% vector of dt's to loop over
dtVec = logspace(smallDt, bigDt, nDt);

% length and width of domain
L = 50;
% base diffusion coefficient
D = 1e0;
% max simulation time
maxT = 10;

%  set D1 and D2, whether isotropic or anisotropic
D1 = D;
D2 = D;
D2 = 0.2 * D;

% define the MT weight functions and analytic solution functions here to
% clean things up within the code
weightFunc = @(X, Y, XX, YY, D, dt, beta)...
               (1 / (beta^-1 * pi * (4 * D * dt))) *...
               exp(-(((XX - X).^2 + (YY - Y).^2) / (beta^-1 * (4 * D * dt))));
weightFuncAniso = @(X, Y, XX, YY, D1, D2, dt, beta)...
                  (1 / (beta^-1 * pi * (4 * sqrt(D1 * D2) * dt))) *...
                  exp(-((XX - X).^2 / (beta^-1 * (4 * D1 * dt)) +...
                  (YY - Y).^2 / (beta^-1 * (4 * D2 * dt))));

analyticFunc = @(X, Y, sigma, D, maxTime)...
               (1 / (2 * pi * (sigma + 2 * D * maxTime))) *...
               exp(-(((0.5 * L - X).^2 + (0.5 * L - Y).^2) /...
               (2 * (sigma + 2 * D * maxTime))));
analyticFuncAniso = @(X, Y, sigma, D1, D2, maxTime)...
                    (1 / (2 * pi * (sigma + 2 * sqrt(D1 * D2) * maxTime))) *...
                    exp(-(((0.5 * L - X).^2 / (2 * (sigma + 2 * D1 * maxTime)) +...
                    (0.5 * L - Y).^2 / (2 * (sigma + 2 * D2 * maxTime)))));

% beta from eq. (6)--controls bandwidth of kernel
beta = 1;

% optimal dt calculation
optDt = (L ./ sqrt(nVec)).^2 ./ ((1 / beta) * 2 * D);
optDtX = (L ./ sqrt(nVec)).^2 ./ ((1 / beta) * 2 * D1);
optDtY = (L ./ sqrt(nVec)).^2 ./ ((1 / beta) * 2 * D2);
optDt = max(optDtX, optDtY);

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
sigma = 1.5e-2 * sqrt(4 * D);

% cell array for holding errors for plotting
errVec = cell(nNum, nDt);

%% run the simulation

% particle number (N) loop
for idxN = 1 : nNum
    
%     print info to screen
    fprintf('N loop: %i of %i\n', idxN, nNum)
   
%     set N
    N = nVec(idxN);
%     space the particles evenly in the x- and y-directions
    Xlin = linspace(0, L, sqrt(N))';
    X = zeros(N, 2);
    [XX, YY] = meshgrid(Xlin, Xlin);
    X(:, 1) = reshape(XX, N, 1);
    X(:, 2) = reshape(YY, N, 1);
    
%     dt loop
    for idxDt = 1 : nDt
        
%         print info to screen
        fprintf('dt loop: %i of %i\n', idxDt, nDt)
        
%         set dt and calculate the number of steps to take
        dt = dtVec(idxDt);
        nSteps = maxT / dt;

%         gaussian IC (assumes domain is [0, L])
        mass = (1 / (2 * pi * sigma^2)) * exp(-((X(:, 1) - 0.5 * L).^2 + (X(:, 2) - 0.5 * L).^2) / (2 * sigma^2));
        initial = mass;
% %         calculate the analytic solution based on istropic/anisotropic
% %         condition
%         analytic = analyticFunc(X(:, 1), X(:, 2), sigma, D, maxT);
        analytic = analyticFuncAniso(X(:, 1), X(:, 2), sigma, D1, D2, maxT);
%         normalize both to sum to 1
        mass = mass / sum(mass);
        analytic = analytic / sum(analytic);

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
% %             calculate the weights based on istropic/anisotropic condition
            val(start : finish) = weightFuncAniso(X(i, 1), X(i, 2), X(idx{i}, 1), X(idx{i}, 2), D1, D2, dt, beta);
%             val(start : finish) = weightFunc(X(i, 1), X(i, 2), X(idx{i}, 1), X(idx{i}, 2), D, dt, beta);
            start = finish + 1;
        end
        
%         when N is super large the sparse methods and clearing large vectors
%         when done saves memory
        clear idx r

%         create the sparse weight matrix
        Wmat = sparse(row, col, val);
        clear row col val
        
%         normalize Wmat such that it is symmetric, using the arithmetic mean of
%         colsum and rowsum--this looks super goofy but is necessary to
%         keep things sparse and save memory
%         Note that it is algebraically equivalent to:
%             Wmat = Wmat ./ ((rowsum + colsum) / 2);
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

%         save the error (and the mass vector, in case you want to plot)
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
if D1 == D2
    strOptDT = '\boldmath $\widehat{\Delta t} = \frac{\left(L / \sqrt{N}\right)^2}{2 D \beta^{-1}}$';
    legend([h fakeDash], [Nlegend ; strOptDT],...
           'location', 'northwest')
else
    strOptDT = '\boldmath $\widehat{\Delta t} = \max\left\{\frac{\left(L / \sqrt{N}\right)^2}{2 D_x \beta^{-1}}, \frac{\left(L / \sqrt{N}\right)^2}{2 D_y \beta^{-1}}\right\}$';
    legend([h fakeDash], [Nlegend ; strOptDT],...
            'location', 'southeast')
end
xlabel('\boldmath$\Delta t$')
ylabel('\textbf{RMSE}')
set(gca,'XDir','reverse','XScale','log','YScale','log');
axis tight