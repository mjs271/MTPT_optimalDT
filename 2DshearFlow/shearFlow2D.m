%==========================================================================
% This script conducts an error analysis of the MTPT method over a range of 
% particle numbers (N) and time step lengths (dt) for a 2D problem with
% linear shear flow.

% This script, will generate Figure 7 in:
%     "Optimal time step length for Lagrangian, interacting-particle simulations
%      of diffusive mixing," XXX Journal XXXX year.
%==========================================================================

%%

clear variables

% weight function for mass transfer
weightFuncAniso = @(X, Y, XX, YY, D1, D2, dt, beta)...
                  (1 / (beta^-1 * pi * (4 * sqrt(D1 * D2) * dt))) *...
                  exp(-((XX - X).^2 / (beta^-1 * (4 * D1 * dt)) +...
                  (YY - Y).^2 / (beta^-1 * (4 * D2 * dt))));


% shear rate in Eq. (20)
al = 0.72;

% max simulation time
t_max=4.2;

% beta from eq. (6)--controls bandwidth of kernel
beta = 1;

% Spacing setup goes here. MUST ensure a single particle resides at
% [x y]=[0 0], then even spacing all around it. Implies np_x,np_y are ODD

% x and y domain limits
x_lims=8*[-1 1]; % MUST be symmetric about the origin, but limits can differ
y_lims=x_lims;

D_XF=0.4;           % Diffusion coefficient for X direction

D_XF2=D_XF;       % Diffusion coefficient for Y direction

% base number of particles
npBase = 51;

% list of particle numbers (discretization in the x- and y-directions)
% Note: this gets added to npBase
np_list=[0 0; 50 0;  0 50; 50 50];

% cell array for holding errors
ErrOut=cell(length(np_list),1);

% number of loops to conduct over N
[np_sz,~]=size(np_list);

% particle number loop
for mm=1:np_sz

%     set number of particles
    np_x=npBase+np_list(mm,1);
    np_y=npBase+np_list(mm,2);
    np=np_x*np_y;

%     check to see if number of particles is valid for this setup
    if mod(np_x,2)==0; error('Need ODD np_x'); end
    if mod(np_y,2)==0; error('Need ODD np_y'); end

%     print some info to screen
    clc;
    fprintf('--> Set number: %i\n',mm);

%     set particle locations
    x_t1=linspace(x_lims(1),0,ceil(np_x/2)); 
    x_t2=linspace(0,x_lims(2),ceil(np_x/2)); 
    x_loc=[x_t1 x_t2(2:end)];   % You get some small round off errors but most all are the same spacing
    dx_xf=mean(diff(x_loc));
    y_t1=linspace(y_lims(1),0,ceil(np_y/2));
    y_t2=linspace(0,y_lims(2),ceil(np_y/2));
    y_loc=[y_t1 y_t2(2:end)];
    dy_xf=mean(diff(y_loc));
%     create a position meshgrid
    [X_loc, Y_loc]=meshgrid(x_loc,y_loc);

%     save initial positions
    Xp0=zeros(np,4);
    Xp0(:,1)=reshape(X_loc,np,1);
    Xp0(:,2)=reshape(Y_loc,np,1);

%     Now shift everybody's initial positions so that the shearing doesn't impact
%     the edges of the particle ensemble as much
    X_loc=X_loc-al.*Y_loc.*t_max/2;
%     put the new positions in the particle array
    Xp=zeros(np,4);
    Xp(:,1)=reshape(X_loc,np,1);
    Xp(:,2)=reshape(Y_loc,np,1);

%     these are percentages of the optimal dt to explore that area
    dt_list=[1.1 1.05 1 0.95 0.9];

%     calculate optimal dt
    optDtX = (16 / np_x)^2 ./ ( 2 * (1 / beta) * 2 * D_XF);
    optDtY = (16 / np_y)^2 ./ ( 2 * (1 / beta) * 2 * D_XF2);
    dtopt = max(optDtX, optDtY);

%     calculate the initial dt's near dtopt, then add some more values to
%     explore the behavior
    DT_list=dtopt.*dt_list;
    DT_list=[DT_list 0.1 0.3 0.4 0.08 0.06 0.05 0.04 0.35 0.03 0.02 0.18 0.125 0.11 0.115 0.105 0.001 0.002];
    DT_list=sort(DT_list,2,'ascend');

%     array to save errors
    Er_mat=zeros(length(dt_list),1);

%     dt loop
    for nn=1:length(DT_list)
%         reset to initial positions
        Xp=zeros(np,3);
        Xp(:,1)=reshape(X_loc,np,1);
        Xp(:,2)=reshape(Y_loc,np,1);

%         set current dt and calculate max time
        dt_xf=DT_list(nn);
        nt_xf=floor(t_max/dt_xf);   % Round nt down a time step if they don't match exactly 
        t_max_xf=nt_xf*dt_xf;

        [Ic,~]=find(Xp(:,1)==0 & Xp(:,2)==0);

        Xp(Ic,3)=1/(dx_xf.*dy_xf);      % Normalize mass by "cell" size

%         search distance for rangesearch
        pad_size=6;
        dzt=pad_size*sqrt(2*(2*D_XF)*dt_xf); 

%         Just focus on the particles within this box to reduce edge impacts
        mon_fac=0.75;   % LESS THAN ONE
        mon_ID=find((Xp0(:,1)>=x_lims(1)*mon_fac)&(Xp0(:,1)<=x_lims(2)*mon_fac)&...
                    (Xp0(:,2)>=y_lims(1)*mon_fac)&(Xp0(:,2)<=y_lims(2)*mon_fac));

%         start timer
        tic
%         set initial simulation time
        t_now=0;

%         time stepping loop
        for nstep=1:nt_xf   

%             forward Euler advection solve
            Xp(:,1)=Xp(:,1) + al.*Xp(:,2).*dt_xf; 

%             conduct the rangesearch to find nearby particles
            [idx, r] = rangesearch(Xp(:,1:2), Xp(:,1:2), 2 * dzt, 'BucketSize',...
                                   ceil(1e-2 * np));

%             determine how many particles are nearby and preallocate the vectors
%             to build the sparse weight matrix
            Nclose = sum(cellfun('length', idx));
            row = zeros(Nclose, 1);
            col = zeros(Nclose, 1);
            val = zeros(Nclose, 1);

%             calculate the entries of the weight matrix
            start = 1;
            for i = 1 : np
                finish = start - 1 + length(idx{i});
%                 evaluate mass-transfer kernel
                row(start : finish) = idx{i};
                col(start : finish) = i;
                val(start : finish) = weightFuncAniso(Xp(i, 1), Xp(i, 2),...
                                      Xp(idx{i}, 1), Xp(idx{i}, 2), D_XF,...
                                      D_XF, dt_xf, beta);
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
            Tmat = speye(np) - spdiags(sum(beta .* Wmat, 2), 0, np, np) + beta .* Wmat;

            clear Wmat

            Xp(:, 3) = Tmat * Xp(:, 3);

            t_now=t_now+dt_xf;

        end
%         stop the timer
        toc

%         Find the analytical solution and the errot (for the inner particles)
        XpM=zeros(length(mon_ID),5);
        XpM(:,1:3)=Xp(mon_ID,1:3);

        for i=1:length(XpM)
            XpM(i,4)=LinShear_full(XpM(i,1),XpM(i,2),D_XF,D_XF2,al,t_now);
        end

        RMSE=sqrt(mean((XpM(:,4)-XpM(:,3)).^2));

        Er_mat(nn)=RMSE;
    end

%     save the errors
    ErrOut{mm}.Er_mat=Er_mat;
    ErrOut{mm}.DT_list=DT_list;
    id=find(DT_list==dtopt);
    ErrOut{mm}.id=id;

end

%% plotting

% set default figure settings and get the color order for plotting
[figSettings, color, rgbcmy, mList, lineList, dashList] = setFigureDefaults();

m_list={'x','d','o','p','+','*','^'};
figure(123);
clf
set(gcf,'color','w');
p_ind=zeros(np_sz,1);
p_names=cell(np_sz,1);
ylms=[0 1e-3];

for nn=1:4

p_ind(nn)=plot(ErrOut{nn}.DT_list,ErrOut{nn}.Er_mat,'-','Marker',m_list{nn});
grid on; hold on;

NPX = np_list(nn,1)+npBase;
NPY = np_list(nn,2)+npBase;

if nn <= 3
    line(dtopt.*[1 1], [((nn - 1) / 3) * ylms(2), (nn / 3) * ylms(2)], 'linestyle', '--','color',get(p_ind(nn),'color'));
else
    line(dtopt.*[1 1],ylms, 'linestyle', '--','color',get(p_ind(nn),'color'));
end
p_names{nn}= ['\boldmath $N_x=$\bf'  num2str(np_list(nn,1)+npBase) ', ' '$N_y=$' num2str(np_list(nn,2)+npBase)];
end
fakeDash = plot(nan, nan, 'k--');
strDTopt = '\boldmath $\widehat{\Delta t} = \max\left\{\frac{\left(L / N_x\right)^2}{4 D \beta^{-1}}, \frac{\left(L / N_y\right)^2}{4 D \beta^{-1}}\right\}$';
legend([p_ind; fakeDash],[p_names; strDTopt],'location','NorthWest');
xlabel('\boldmath $\Delta t$'); ylabel('\textbf{RMSE}');
set(gca,'XDir','reverse','XScale','log');
ylim(ylms);
