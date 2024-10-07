
clear all

data_list = [{'img_sos_tr6_n5'},...
    {'img_sos_tr6_n6'},...
    {'img_sos_tr6_n7'},...
    {'img_sos_tr6_n12'},...
    {'img_sos_me_flash'}];

for data_idx=1:length(data_list)
    load(strcat("../Data/",data_list{data_idx},".mat"))
    
    n_echoes = 5; % How many echoes to fit?
    TE_list=3:6:((6*n_echoes)-3); % Effective echo times in n-periodic sequences
    
    % Initial guess
    T2s0=15.0; 
    S0=1.5; 
    lb=[0 0];
    ub=[Inf 500];
    ms_pts = 50;    % 50 multistart pts - reduce this if you want to fit more quickly!
    
    
    S0fit=zeros(256,256);
    T2sfit=zeros(256,256);
    
    mask=(img_sos(:,:,1)>(0.05*max(img_sos(:)))); % 5 percent of max for mask
    %%
    tic 
    parpool;
    for row=1:256
        for col=1:256
            
                if mask(row,col)>0
                    t2s_fitfun = @(x,TE_list)x(1)*exp(-TE_list/x(2));
                    x0=[img_sos(row,col,1)*S0*1e6,T2s0];
    
                    problem = createOptimProblem('lsqcurvefit','x0',x0, ...
                        'objective',t2s_fitfun,'lb',lb,'ub',ub, ...
                        'xdata', TE_list(1:n_echoes), ...
                        'ydata',squeeze(img_sos(row,col,1:n_echoes)).'*1e6);
    
                    ms = MultiStart('UseParallel',true,'Display',"off");
    
                    [xfit,~] = run(ms,problem,ms_pts);
    
                    T2sfit(row,col)=xfit(2);
                    S0fit(row,col)=xfit(1);
                end       
            
            
        end

        imagesc(rot90(mask./T2sfit,2));
        axis image off
        clim([0 0.15])
        colormap parula
        drawnow
    
    end
    delete(gcp)
    toc
    
    save(strcat("../Data/FitPars_",data_list{data_idx},".mat"),"T2sfit","S0fit");
end