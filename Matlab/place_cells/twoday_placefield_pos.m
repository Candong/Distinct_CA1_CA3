clear all; %close all;

%%
[PC_filepaths, temp]=uigetfile('*.mat', 'Chose f & n files to load:','MultiSelect','on');
planenum = input(['What are the planes you would like to analysis  >> ']);
for f=1:size(PC_filepaths,2)
 if contains(PC_filepaths{f},'_f_')

    for p=planenum
    if contains(PC_filepaths{f},['plain' num2str(p)]) 
        load(PC_filepaths{f});        
        f_sig_PFs{p}=sig_PFs;
        f_sig_PFs_with_noise{p}=sig_PFs_with_noise;
        f_meantrans_withnoise{p}=mean_trans;
        f_PF_start_bins{p}=PF_start_bins;
        f_PF_end_bins{p}=PF_end_bins;
            
    end
    end
 elseif contains(PC_filepaths{f},'_n_')

    for p=planenum
    if contains(PC_filepaths{f},['plain' num2str(p)]) 
        load(PC_filepaths{f});
        
        n_sig_PFs{p}=sig_PFs;
        n_sig_PFs_with_noise{p}=sig_PFs_with_noise;
        n_meantrans_withnoise{p}=mean_trans;
        n_PF_start_bins{p}=PF_start_bins;
        n_PF_end_bins{p}=PF_end_bins;
            
    end  
  end  
   
 end
end

%% 
total_cor=[];
exclude_total_cor=[];
common_total_cor=[];
fn_total_cor=[];
for p=planenum
  
%   for i=1:size(f_sig_PFs{p},2)
[f_pf_num{p} f_pf_id{p}]=find(cellfun(@isempty,f_sig_PFs{p})==0);
[n_pf_num{p} n_pf_id{p}]=find(cellfun(@isempty,n_sig_PFs{p})==0);
      
f_n_PCid{p}=union(f_pf_id{p},n_pf_id{p});
cur_f_mean=f_meantrans_withnoise{p};
cur_n_mean=n_meantrans_withnoise{p};
transient_cor{p}=diag(corr(cur_f_mean(f_n_PCid{p},:)',cur_n_mean(f_n_PCid{p},:)'));
      
total_cor=[total_cor transient_cor{p}'];   

exclude_PC_id{p}=setxor(f_pf_id{p},n_pf_id{p});
exclude_cor{p}=diag(corr(cur_f_mean(exclude_PC_id{p},:)',cur_n_mean(exclude_PC_id{p},:)'));
exclude_total_cor=[exclude_total_cor exclude_cor{p}'];    


common_PC_id{p}=intersect(f_pf_id{p},n_pf_id{p});
if ~isempty(common_PC_id{p})
common_cor{p}=diag(corr(cur_f_mean(common_PC_id{p},:)',cur_n_mean(common_PC_id{p},:)'));
common_total_cor=[common_total_cor common_cor{p}']; 
end
 
common_PC_id{p}=intersect(f_pf_id{p},n_pf_id{p});
if ~isempty(common_PC_id{p})
common_cor{p}=diag(corr(cur_f_mean(common_PC_id{p},:)',cur_n_mean(common_PC_id{p},:)'));
common_total_cor=[common_total_cor common_cor{p}']; 
end


if ~isempty(f_pf_id{p})
fn_cor{p}=diag(corr(cur_f_mean(f_pf_id{p},:)',cur_n_mean(f_pf_id{p},:)'));
fn_total_cor=[fn_total_cor fn_cor{p}']; 
end
    
end
figure;
plot([1],total_cor,'o')
hold on;plot([1],mean(total_cor),'*')
figure;
plot([1],exclude_total_cor,'o')
hold on;plot([1],mean(exclude_total_cor),'*')


figure;
plot([1],common_total_cor,'o')
hold on;plot([1],mean(common_total_cor),'*')


figure;
plot([1],fn_total_cor,'o')
hold on;plot([1],mean(fn_total_cor),'*')

range=-1:0.1:1;
[fn_dist edges]=histcounts(fn_total_cor,range);
figure;
histogram(fn_total_cor,range)
ylim([0 10])

%%  


% [meanCOM cur_SP]=meanCOMandSP(PF_denoise,start_b,end_b,numbins);


%%
function [meanCOM SP]=meanCOMandSP(binmean,start_b,end_b,numbins)
binM=1:numbins;
COM=sum(binmean(:,start_b:end_b).*binM(start_b:end_b),2)./sum(binmean(:,start_b:end_b),2);
COM(isnan(COM))=0;
A=max(binmean(:,start_b:end_b),[],2);
meanCOM=sum(A.*COM)/sum(A);

SP=1/(sqrt(sum(A.*((COM-meanCOM).^2))/sum(A)));
 end 