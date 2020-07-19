clear all; %close all;
[behavior, temp]=uigetfile('*.mat', 'Chose f and n files to load:','MultiSelect','on');
file_count=1;

if ~isa(behavior,'cell') 
  behavior={behavior};  
    
end 
for f=1:size(behavior,2)
    if contains(behavior{f},'_f_')
       behavior_filepaths{file_count}= behavior{f};
       file_count=file_count+1;
    
    end
end

behavior_filepaths=sort(behavior_filepaths);

%%
f_total_act=[];
f_total_infield_act=[];
f_total_outfield_act=[];
for f=1:size(behavior_filepaths,2)
load([temp behavior_filepaths{f}]);

pf_id=find(~isnan(number_of_PFs));

for i=1:length(pf_id)  
    binmean_temp=sig_PFs_with_noise{1,pf_id(i)};
    total_bins=size(binmean_temp,1)*size(binmean_temp,2);
    total_act=sum(sum(binmean_temp));
    infield_act=0;
   
    for j=1:number_of_PFs(pf_id(i))
            
            binmean_PF_temp=sig_PFs{j,pf_id(i)};            
            mean_PF=mean(binmean_PF_temp,2);
            mean_PF_withnoise=mean(binmean_temp,2); 
            start_bin=PF_start_bins(j,pf_id(i));
            end_bin=PF_end_bins(j,pf_id(i));
            infield_act=infield_act+sum(sum(binmean_temp(start_bin:end_bin,:)));
            infield_bins=(end_bin-start_bin+1)*size(binmean_temp,2);
    end
    outfield_act=total_act-infield_act;
    f_total_act=[f_total_act total_act/total_bins];
    f_total_infield_act=[f_total_infield_act infield_act/(infield_bins)];   
    f_total_outfield_act=[f_total_outfield_act outfield_act/(total_bins-infield_bins)];
end


end

 %% 
% %clear all
% 
% 
%[behavior, temp]=uigetfile('*.mat', 'Chose f and n files to load:','MultiSelect','on');
file_count=1;

if ~isa(behavior,'cell') 
  behavior={behavior};  
    
end 
for f=1:size(behavior,2)
    if contains(behavior{f},'_n_')
       behavior_filepaths{file_count}= behavior{f};
       file_count=file_count+1;
    
    end
end

behavior_filepaths=sort(behavior_filepaths);

n_total_act=[];
n_total_infield_act=[];
n_total_outfield_act=[];
for f=1:size(behavior_filepaths,2)
load([temp behavior_filepaths{f}]);

pf_id=find(~isnan(number_of_PFs));

for i=1:length(pf_id)  
    binmean_temp=sig_PFs_with_noise{1,pf_id(i)};
    total_bins=size(binmean_temp,1)*size(binmean_temp,2);
    total_act=sum(sum(binmean_temp));
    infield_act=0;
   
    for j=1:number_of_PFs(pf_id(i))
            
            binmean_PF_temp=sig_PFs{j,pf_id(i)};            
            mean_PF=mean(binmean_PF_temp,2);
            mean_PF_withnoise=mean(binmean_temp,2); 
            start_bin=PF_start_bins(j,pf_id(i));
            end_bin=PF_end_bins(j,pf_id(i));
            infield_act=infield_act+sum(sum(binmean_temp(start_bin:end_bin,:)));
            infield_bins=(end_bin-start_bin+1)*size(binmean_temp,2);
    end
    outfield_act=total_act-infield_act;
    n_total_act=[n_total_act total_act/total_bins];
    n_total_infield_act=[n_total_infield_act infield_act/(infield_bins)];   
    n_total_outfield_act=[n_total_outfield_act outfield_act/(total_bins-infield_bins)];
end


end

%%
f_ratio=f_total_outfield_act./f_total_infield_act;
n_ratio=n_total_outfield_act./n_total_infield_act;
g1 = repmat({'F'},length(f_total_infield_act),1);
g2 = repmat({'N'},length(n_total_infield_act),1);
%g3 = repmat({'FN'},length(total_cor),1);
g=[g1; g2];
figure
boxplot([f_ratio,n_ratio],g)

%%
save('CA1_inout_ratio_fn','f_ratio','n_ratio')