%% STEP 1: Reconstructing individual echo images from raw data
% Here we take the raw .dat files (provided separately - download these
% and move them into the 'data' folder) and:
% i) Do a sum-of-squares recon of the FLASH data across coils
% ii) Do a POCS Partial Fourier recon of the individual SSFP echoes for
% each coil
% iii) Combine the per-echo data across coils via sum-of-squares


clear all;

%% Add external functions
startdir=pwd;
addpath(genpath("..\external"))
cd('../data')
%% Read in Multiecho FLASH data
%read in data    
twix = mapVBVD(53);
d_orig=twix.image(:,:,:,:);

%reorder data axes
d_srt=permute(d_orig,[1 3 2 4]);

%correct for raw data oversampling in FE direction
d_srt = d_srt(1:2:end,:,:,:);

%% Compute sum-of-squares across coils for each echo and write mat files
%compute sum-of-squares across coils
for echo=1:size(d_srt,4)
    img_sos(:,:,echo)=fliplr(sos(ifft2c(d_srt(:,:,:,echo))));
end

% Save .mat file
save(strcat('img_sos_me_flash.mat'),"img_sos");

cd(startdir)

%% Loop through 4 SSFP datasets
for cont_idx_=1:4

    clearvars -except startdir cont_idx_
    if cont_idx_==1
        % file no 57, N=12, TR=6ms, grad spoiling = 20 lines, central
        % config state = F5
        cont_list = [57]; N=12; TR=6; spl_lines=20; cfg=5;
    elseif cont_idx_==2
        % file no 58, N=6, TR=6ms, grad spoiling = 40 lines, central
        % config state = F2        
        cont_list = [58]; N=6; TR=6; spl_lines=40; cfg=2;
    elseif cont_idx_==3
        % file no 59, N=5, TR=6ms, grad spoiling = 40 lines, central
        % config state = F2        
        cont_list = [59]; N=5; TR=6; spl_lines=40; cfg=2;
    elseif cont_idx_==4
        % file no 60, N=7, TR=6ms, grad spoiling = 34 lines, central
        % config state = F3       
        cont_list = [60]; N=7; TR=6; spl_lines=34; cfg=3;
    end
        
    %read in data        
    ncont=length(cont_list);
    for cont_no = 1:length(cont_list)
        twix = mapVBVD(cont_list(cont_no));
        d_orig=twix.image.unsorted;
    end

    %reorder data axes
    for echo=1:N
        d_pss_n(:,:,:,echo) = permute(d_orig(1:2:end,:,echo:N:end),[1 3 2]);
    end
    % How many coils were used?
    ncoil=32;

    %% Compute each of the F-state signals (demodulate by phase term)
    fstates=[0:N-1];
    temp_idx=0;
    for f_idx=fstates
        temp_idx=temp_idx+1;
        d_fstates_n(:,:,:,temp_idx) = mean(d_pss_n.*reshape(exp(f_idx*1i*[0:(N-1)]/N*2*pi),1,1,1,[]),4);
    end

    %% Perform partial Fourier reconstructions with POCS
    % see https://ece-classes.usc.edu/ee591/library/Pauly-PartialKspace.pdf

    iters=10; % choose number of POCS iterations
    slc=1;    % choose slice (we only have one slice)

    %Loop through each of the N F-states
    for f_idx=1:N

        %Create empty array (nFE,nPE,nCoils) for POCS results
        S = zeros(256,256,ncoil);

        %Work out which lines we sample, and which we don't (which POCS will
        %calculate)
        pe_mask = circshift([ones(1,256),zeros(1,256)],spl_lines*(f_idx-cfg-1));
        pe_lines = find(pe_mask(1:256));
        non_pe_lines = setxor(pe_lines,1:256);

        %Choose the central lines around each echo to work out the image phase
        mid_mask = zeros(size(S,[1 2]));
        sym_sz = sum(pe_mask(1:256))-128-1;
        mid_mask(:,128-sym_sz:128+sym_sz-1)=1;
        n_lines = length(pe_lines);

        %Fill in measured lines into S
        if f_idx>(cfg+1)  % Is F-state to left or right of centre?
            S(:,pe_lines,:)=d_fstates_n(:,1:n_lines,:,f_idx);    
        else 
            S(:,pe_lines,:)=d_fstates_n(:,end-n_lines+1:end,:,f_idx);    
        end

        % Compute phase consistency term for POCS
        phase_im = exp(1i.*angle(ifft2c(mid_mask.*S)));

        %% Run POCS reconstruction for each coil and visualise results
        figure(1);clf   % clear figure for iteration results

        for idx_c=1:ncoil

            % fill k-space with known data
            k_pocs = zeros(size(S,[1 2]));             
            k_pocs(:,pe_lines)=S(:,pe_lines,idx_c); 

            for it = 1:iters   %for each POCS iteration

                im_pocs = abs(ifft2c(k_pocs)).*phase_im(:,:,idx_c);    % enforce phase consistency
                k_new = fft2c(im_pocs);                                % update k-space guess
                k_pocs(:,non_pe_lines)=k_new(:,non_pe_lines);          % combine with known data


                % show results during POCS iterations
                figure(1)
                subplot(2,2,1) % Zerofilled k-space: symmetrically sampled region
                imagesc(log(abs(mid_mask.*k_pocs))); axis off equal
                title("ZF - slice "+num2str(slc)+" coil "+num2str(idx_c))
                CX = caxis;

                subplot(2,2,2) % Zerofilled image: symmetrically sampled region
                imagesc(abs(ifft2c(mid_mask.*k_pocs))); axis off square
                title("ZF - slice "+num2str(slc)+" coil: "+num2str(idx_c))

                subplot(2,2,3) % Reconstructed data: completed with POCS
                imagesc(log(abs(k_pocs))); axis off equal
                title("SR - slice "+num2str(slc)+" coil "+num2str(idx_c)+" iter: "+num2str(it))
                caxis(CX);

                subplot(2,2,4) % Reconstructed image: completed with POCS
                imagesc(abs(im_pocs)); axis off square
                title("SR - slice "+num2str(slc)+" coil "+num2str(idx_c)+" iter: "+num2str(it))

                colormap gray
                drawnow

            end %ending iteration loop for POCS

            % store final k-space results per coil
            k_fin(:,:,idx_c,f_idx)=single(k_pocs);

        end  %ending coil loop for POCS
    end %ending F-state loop for this dataset

    %% Compute sum-of-squares across coils for each F-state and write DICOM files
    clear img_sos
    %compute sum-of-squares across coils
    for echo=1:N % #Echo=#Fstate+1   (i.e. echo1 = F0 signal)
        img_sos(:,:,echo)=fliplr(sos(ifft2c(k_fin(:,:,:,echo))));
    end

    save(strcat('img_sos_tr',num2str(TR),'_n',num2str(N),'.mat'),"img_sos");

end
