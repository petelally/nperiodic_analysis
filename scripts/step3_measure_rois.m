

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

%%
load('T2sS0_fit_me_tr6_fa90.mat')
R2s_gt_rois = rois.*fliplr(1./T2sfit);
R2s_gt_means = nansum(reshape(R2s_gt_rois,256*256,[]),1)./sum(roi_mask(:));
R2s_gt_std = nanstd(reshape(R2s_gt_rois,256*256,[]),0,1);

figure(1)
clf 
subplot(1,5,1)
imagesc(img_sos(:,:,4));axis off image;colormap gray;
figure(2)
subplot(1,5,1)
imagesc(1./T2sfit.*(T2sfit>0));axis off image;colormap parula;clim([0 0.15])

%%
load('../Data/FitPars_img_sos_tr6_n5.mat')
R2s_n5_rois = rois.*fliplr(1./T2sfit);
R2s_n5_means = nansum(reshape(R2s_n5_rois,256*256,[]),1)./sum(roi_mask(:));
R2s_n5_std = nanstd(reshape(R2s_n5_rois,256*256,[]),0,1);

figure(1)
subplot(1,5,2)
imagesc(img_sos(:,:,4));axis off image;colormap gray;
figure(2)
subplot(1,5,2)
imagesc(1./T2sfit.*(T2sfit>0));axis off image;colormap parula;clim([0 0.15])

figure(3)
clf
errorbar(squeeze(R2s_gt_means),(squeeze(R2s_gt_means)-squeeze(R2s_n5_means))./squeeze(R2s_gt_means)*100,...
    squeeze(R2s_n5_std)./squeeze(R2s_gt_means)*100,'or','MarkerFaceColor','r')
%%
load('T2sS0_fit_tr6_n6.mat')
R2s_n6_rois = rois.*fliplr(1./T2sfit);
R2s_n6_means = nansum(reshape(R2s_n6_rois,256*256,[]),1)./sum(roi_mask(:));
R2s_n6_std = nanstd(reshape(R2s_n6_rois,256*256,[]),0,1);

figure(1)
subplot(1,5,3)
imagesc(img_sos(:,:,4));axis off image;colormap gray;
figure(2)
subplot(1,5,3)
imagesc(1./T2sfit.*(T2sfit>0));axis off image;colormap parula;clim([0 0.15])

figure(3)
hold on
errorbar(squeeze(R2s_gt_means)+0.001,(squeeze(R2s_gt_means)-squeeze(R2s_n6_means))./squeeze(R2s_gt_means)*100,...
    squeeze(R2s_n6_std)./squeeze(R2s_gt_means)*100,'ob','MarkerFaceColor','b')
%%
load('T2sS0_fit_tr6_n7.mat')
R2s_n7_rois = rois.*fliplr(1./T2sfit);
R2s_n7_means = nansum(reshape(R2s_n7_rois,256*256,[]),1)./sum(roi_mask(:));
R2s_n7_std = nanstd(reshape(R2s_n7_rois,256*256,[]),0,1);

figure(1)
subplot(1,5,4)
imagesc(img_sos(:,:,4));axis off image;colormap gray;
figure(2)
subplot(1,5,4)
imagesc(1./T2sfit.*(T2sfit>0));axis off image;colormap parula;clim([0 0.15])

figure(3)
hold on
errorbar(squeeze(R2s_gt_means)+0.002,(squeeze(R2s_gt_means)-squeeze(R2s_n7_means))./squeeze(R2s_gt_means)*100,...
    squeeze(R2s_n7_std)./squeeze(R2s_gt_means)*100,'om','MarkerFaceColor','m')
%%
load('T2sS0_fit_tr6_n12.mat')
R2s_n12_rois = rois.*fliplr(1./T2sfit);
R2s_n12_means = nansum(reshape(R2s_n12_rois,256*256,[]),1)./sum(roi_mask(:));
R2s_n12_std = nanstd(reshape(R2s_n12_rois,256*256,[]),0,1);

figure(1)
subplot(1,5,5)
imagesc(img_sos(:,:,4));axis off image;colormap gray;set(gcf,'Color','w')
figure(2)
subplot(1,5,5)
imagesc(1./T2sfit.*(T2sfit>0));axis off image;colormap parula;clim([0 0.15]);set(gcf,'Color','w')

figure(3)
hold on
errorbar(squeeze(R2s_gt_means)+0.003,(squeeze(R2s_gt_means)-squeeze(R2s_n12_means))./squeeze(R2s_gt_means)*100,...
    squeeze(R2s_n12_std)./squeeze(R2s_gt_means)*100,'ok','MarkerFaceColor','k')
%%
plot([0 0.2],[0 0],'k:')
set(gcf,'Color','w')
ylabel("Estimated $R_{2}^*$ vs Multi-echo FLASH $R_{2}^*$ ($\%$)",Interpreter='latex')
xlabel("Multi-echo FLASH $R_{2}^*$ ($ms^{-1}$)",Interpreter='latex')
set(gca,'TickLabelInterpreter','latex','FontSize',14);

