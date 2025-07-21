
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Rahul Gaurav and Romain Valabregue
% Date: 2016
% Last Modified: 21st July 2025
% Description: This script performs analysis of substantia nigra-based normalized signal intensity (NSI).
% Version: 2.0
% Contact: rahul.gaurav@icm-instutute.org
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function SN_NM_NSI(suj, pool, TheFile, fvol, fref, froi_total)
progressbar('Nigral Normalized Signal Intensity Analysis')

cout = struct;
[pp sujn] = get_parent_path(suj);
cout.pool = pool
cout.suj = sujn
for k=1:length(fvol)
    
    progressbar(k/length(fvol)) % Update progress bar
    [h Aref] = nifti_spm_vol(fref{k});
    
    [g As]  = nifti_spm_vol(fvol{k});
    volSr1L=0; volSr1R=0; % Trial: RG
    [h Ar1L] = nifti_spm_vol(froi_total{k});
    
    non_zero_ind_ref = Aref>0;
    non_zero_ind_roiL = Ar1L>0;
    
    %take the left and the right roi and find the x mm position at the middle of both
    % from index to mm
    ind = [ir jr kr ones(size(ir))]';     indxyzL = h.mat*ind;
   
    RLseparation = mean(indxyzL(1,:)]);
    
    [ir, jr, kr] = ind2sub(size(Aref),find(Aref>0));
    ind = [ir jr kr ones(size(ir))]'; indxyz = h.mat*ind;
    indleft = indxyz(1,:)<= RLseparation;
    
    indmatL = sub2ind(size(Aref),ir(indleft),jr(indleft),kr(indleft));
  
    
    %Slice by slice
    
    nslice = unique(kr); %z index of the reference region
    SratioR=[];SratioL=[]; volSr1R=[];volSr1L=[];
    S_CNR_R=[];S_CNR_L=[];
    
    for nbs = 1:length(nslice) 
        As_slice = As(:,:,nslice(nbs));
        Aref_slice = Aref(:,:,nslice(nbs));

        non_zero_ind_ref  = Aref_slice>0;

        Sbnd_bg_mean = mean(As_slice(non_zero_ind_ref));
        Sbnd_bg_std  = std (As_slice(non_zero_ind_ref));
        
        Sr1R =  mean(As_slice(non_zero_ind_roiR));     Sr1L =  mean(As_slice(non_zero_ind_roiL));
        volSr1R(end+1) = length(As_slice(non_zero_ind_roiR)); volSr1L(end+1) = length(As_slice(non_zero_ind_roiL));
        SratioR(end+1) = Sr1R ./Sbnd_bg_mean*100; SratioL(end+1) = Sr1L ./Sbnd_bg_mean*100;
        S_CNR_R(end+1) = (Sr1R - Sbnd_bg_mean)./Sbnd_bg_std ;  S_CNR_L(end+1) = (Sr1L - Sbnd_bg_mean)./Sbnd_bg_std ;
    end
    
    cout.sig_SNL_sl_mean(k) = sum(SratioL.*volSr1L/sum(volSr1L)); % mean(SratioL);
    cout.sig_SNR_sl_mean(k)      = sum(SratioR.*volSr1R/sum(volSr1R)); %mean(SratioR);
    
    cout.sig_CNR_L(k)       = sum(S_CNR_L.*volSr1L/sum(volSr1L)); % mean(SratioL);
    cout.sig_CNR_R(k)            = sum(S_CNR_R.*volSr1R/sum(volSr1R)); %mean(SratioR);
    
    cout.sig_SN_both(k) = (cout.sig_SNL_sl_mean(k) +  cout.sig_SNR_sl_mean(k)) /2 ;
    cout.sig_CNR_both(k) = (cout.sig_CNR_L(k) +  cout.sig_CNR_R(k)) /2 ;
    
    disp(sujn(k))
    disp(k)
    
end

write_result_to_csv(cout,TheFile)
end