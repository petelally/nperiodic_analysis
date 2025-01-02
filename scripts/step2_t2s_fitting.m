%% STEP 2: Generating T2* estimates from sum-of-squares data
% Here we take the sum-of-squares images from Step 1 and fit them to a
% monoexponential decay function to estimate T2* and intercept S0
%
% N.B. Fitting takes ~3h on a quad-core system with 50x multistart points.
% If you want to speed this up, you can reduce the number of multistart
% points (or just do a single fit), and parallelise the fitting across
% different voxels

%% Load in each dataset in a loop for fitting
clear all

data_list = [{'img_sos_tr6_n5'},...
    {'img_sos_tr6_n6'},...
    {'img_sos_tr6_n7'},...
    {'img_sos_tr6_n12'},...
    {'img_sos_me_flash'}];

for data_idx=1:length(data_list)
    load(strcat("../Data/",data_list{data_idx},".mat")); % Load dataset
    
    n_echoes = 5; % How many echoes to fit? We choose 5 echoes
    TE_list=3:6:((6*n_echoes)-3); % Effective echo times for all approaches
    
    % Initial guess
    T2s0=15.0; % 15ms starting guess
    S0=1.5; % Starting assumption intercept is 1.5x signal at TE1
    lb=[0 0]; % Non-negative
    ub=[Inf 500]; % Max T2 is 500 ms
    ms_pts = 50;    % 50 multistart pts - reduce this if you want to fit more quickly!
    
    %Initialise parameter arrays
    S0fit=zeros(256,256);
    T2sfit=zeros(256,256);
    
    %Create mask
    mask=(img_sos(:,:,1)>(0.05*max(img_sos(:)))); % 5 percent of max for mask
    %% Do the fitting!
    tic 
    parpool;
    for row=1:256
        for col=1:256
            
                if mask(row,col)>0
                    t2s_fitfun = @(x,TE_list)x(1)*exp(-TE_list/x(2)); % Monoexponential decay fit
                    x0=[img_sos(row,col,1)*S0*1e6,T2s0]; % Initial guess (rescale magnitude by 1e6 to consistent parameter scales)
    
                    problem = createOptimProblem('lsqcurvefit','x0',x0, ...
                        'objective',t2s_fitfun,'lb',lb,'ub',ub, ...
                        'xdata', TE_list(1:n_echoes), ...
                        'ydata',squeeze(img_sos(row,col,1:n_echoes)).'*1e6);
    
                    ms = MultiStart('UseParallel',true,'Display',"off");
    
                    [xfit,~] = run(ms,problem,ms_pts); % Fit voxel data to decay fn
    
                    T2sfit(row,col)=xfit(2); % Store T2* estimate
                    S0fit(row,col)=xfit(1); % Store S0 estimate
                end       
            
            
        end

        imagesc(rot90(mask./T2sfit,2)); % Display R2* map as it updates row-by-row in real time 
        axis image off
        clim([0 0.15])  % Uncomment if using MATLAB R2022a or later
        % caxis([0 0.15])  % Uncomment if using MATLAB versions prior to R2022a
        colormap parula
        drawnow
    
    end
    delete(gcp)
    toc
    
    save(strcat("../Data/FitPars_",data_list{data_idx},".mat"),"T2sfit","S0fit"); % Save fitted parameter maps
end
