%% STEP 3: Generating ROI statistics
% Here we take the fitted datasets and apply ROI statistics to replucate
% Figure 4 from the paper

%% Create ROIs
% Create a circular mask for the ROI
roi_mask = zeros(256,256);
for x_idx = 1:256
    for y_idx = 1:256
        if(sqrt((x_idx-128).^2+(y_idx-128).^2)<7)
            roi_mask(y_idx,x_idx) = 1;
        end
    end
end


% Positions of the NIST spheres in the phantom images
roi_locs = [-48, 2 ;...
            -39, 32;...
            -14, 49;...
             16, 50;...
             41, 31;...
             52,  3;...
             44,-26;...
             18,-46;...
            -13,-46;...
            -38,-27;...
            -18,-19;...
            -19, 21;...
             22, 23;...
             22,-17];

% Shift the circular ROI to each position 
for roi_no = 1:size(roi_locs,1)
    rois(:,:,roi_no) = circshift(roi_mask,roi_locs(roi_no,:));
end

%% Display and measure the fitted FLASH data ('ground truth')
load('../data/FitPars_img_sos_me_flash.mat')
R2s_gt_rois = rois.*fliplr(1./T2sfit);
R2s_gt_means = nansum(reshape(R2s_gt_rois,256*256,[]),1)./sum(roi_mask(:));
R2s_gt_std = nanstd(reshape(R2s_gt_rois,256*256,[]),0,1);
figure(1)
clf 
subplot(1,5,1)
imagesc(img_sos(:,:,4));axis off image;colormap gray;
title('T2*w image (FLASH)')
figure(2)
subplot(1,5,1)
imagesc(1./T2sfit.*(T2sfit>0));axis off image;colormap parula;
clim([0 0.15]);  % Uncomment if using MATLAB R2022a or later
% caxis([0 0.15]);  % Uncomment if using MATLAB versions prior to R2022a
title('R2* map (FLASH)')

%% Display and measure the fitted SSFP (N=5) data
load('../Data/FitPars_img_sos_tr6_n5.mat')
R2s_n5_rois = rois.*fliplr(1./T2sfit);
R2s_n5_means = nansum(reshape(R2s_n5_rois,256*256,[]),1)./sum(roi_mask(:));
R2s_n5_std = nanstd(reshape(R2s_n5_rois,256*256,[]),0,1);
figure(1)
subplot(1,5,2)
imagesc(img_sos(:,:,4));axis off image;colormap gray;
title('T2*w image (SSFP, N=5)')
figure(2)
subplot(1,5,2)
imagesc(1./T2sfit.*(T2sfit>0));axis off image;colormap parula;
clim([0 0.15]);  % Uncomment if using MATLAB R2022a or later
% caxis([0 0.15]);  % Uncomment if using MATLAB versions prior to R2022a
title('R2* map (SSFP, N=5)')
figure(3)
clf
errorbar(squeeze(R2s_gt_means),(squeeze(R2s_gt_means)-squeeze(R2s_n5_means))./squeeze(R2s_gt_means)*100,...
    squeeze(R2s_n5_std)./squeeze(R2s_gt_means)*100,'or','MarkerFaceColor','r')

%% Display and measure the fitted SSFP (N=6) data
load('../data/FitPars_img_sos_tr6_n6.mat')
R2s_n6_rois = rois.*fliplr(1./T2sfit);
R2s_n6_means = nansum(reshape(R2s_n6_rois,256*256,[]),1)./sum(roi_mask(:));
R2s_n6_std = nanstd(reshape(R2s_n6_rois,256*256,[]),0,1);
figure(1)
subplot(1,5,3)
imagesc(img_sos(:,:,4));axis off image;colormap gray;
title('T2*w image (SSFP, N=6)')
figure(2)
subplot(1,5,3)
imagesc(1./T2sfit.*(T2sfit>0));axis off image;colormap parula;
clim([0 0.15]);  % Uncomment if using MATLAB R2022a or later
% caxis([0 0.15]);  % Uncomment if using MATLAB versions prior to R2022a
title('R2* map (SSFP, N=6)')
figure(3)
hold on
errorbar(squeeze(R2s_gt_means)+0.001,(squeeze(R2s_gt_means)-squeeze(R2s_n6_means))./squeeze(R2s_gt_means)*100,...
    squeeze(R2s_n6_std)./squeeze(R2s_gt_means)*100,'ob','MarkerFaceColor','b')

%% Display and measure the fitted SSFP (N=7) data
load('../data/FitPars_img_sos_tr6_n7.mat')
R2s_n7_rois = rois.*fliplr(1./T2sfit);
R2s_n7_means = nansum(reshape(R2s_n7_rois,256*256,[]),1)./sum(roi_mask(:));
R2s_n7_std = nanstd(reshape(R2s_n7_rois,256*256,[]),0,1);
figure(1)
subplot(1,5,4)
imagesc(img_sos(:,:,4));axis off image;colormap gray;
title('T2*w image (SSFP, N=7)')
figure(2)
subplot(1,5,4)
imagesc(1./T2sfit.*(T2sfit>0));axis off image;colormap parula;
clim([0 0.15]);  % Uncomment if using MATLAB R2022a or later
% caxis([0 0.15]);  % Uncomment if using MATLAB versions prior to R2022a
title('R2* map (SSFP, N=7)')
figure(3)
hold on
errorbar(squeeze(R2s_gt_means)+0.002,(squeeze(R2s_gt_means)-squeeze(R2s_n7_means))./squeeze(R2s_gt_means)*100,...
    squeeze(R2s_n7_std)./squeeze(R2s_gt_means)*100,'om','MarkerFaceColor','m')

%% Display and measure the fitted SSFP (N=12) data
load('../data/FitPars_img_sos_tr6_n12.mat')
R2s_n12_rois = rois.*fliplr(1./T2sfit);
R2s_n12_means = nansum(reshape(R2s_n12_rois,256*256,[]),1)./sum(roi_mask(:));
R2s_n12_std = nanstd(reshape(R2s_n12_rois,256*256,[]),0,1);
figure(1)
subplot(1,5,5)
imagesc(img_sos(:,:,4));axis off image;colormap gray;set(gcf,'Color','w')
title('T2*w image (SSFP, N=12)')
figure(2)
subplot(1,5,5)
imagesc(1./T2sfit.*(T2sfit>0));axis off image;colormap parula;
clim([0 0.15]);  % Uncomment if using MATLAB R2022a or later
% caxis([0 0.15]);  % Uncomment if using MATLAB versions prior to R2022a
set(gcf,'Color','w')
title('R2* map (SSFP, N=12)')
figure(3)
hold on
errorbar(squeeze(R2s_gt_means)+0.003,(squeeze(R2s_gt_means)-squeeze(R2s_n12_means))./squeeze(R2s_gt_means)*100,...
    squeeze(R2s_n12_std)./squeeze(R2s_gt_means)*100,'ok','MarkerFaceColor','k')

%% Annotate the display of the Bland-Altman plot
plot([0 0.2],[0 0],'k:')
set(gcf,'Color','w')
ylabel("Estimated $R_{2}^*$ vs Multi-echo FLASH $R_{2}^*$ ($\%$)",Interpreter='latex')
xlabel("Multi-echo FLASH $R_{2}^*$ ($ms^{-1}$)",Interpreter='latex')
set(gca,'TickLabelInterpreter','latex','FontSize',14);

