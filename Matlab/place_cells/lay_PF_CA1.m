clear all; %close all;
[behavior_filepaths, temp]=uigetfile('*.mat', 'Chose f files to load:','MultiSelect','on');
behavior_filepaths=sort(behavior_filepaths);
PF_center=[];
PF_mean=[];
PF_count=1;
placecell_id={};
MEAN=[];
for f=1:size(behavior_filepaths,2)
load([temp behavior_filepaths{f}]);
placecell_id{f}=[];
for i=1:size(sig_PFs,2)
    for j=1:size(sig_PFs,1)
        if ~isempty(sig_PFs{j,i})
            binmean_temp=sig_PFs_with_noise{j,i};
            binmean_temp(isnan(binmean_temp))=0;
            binmean_PF_temp=sig_PFs{j,i};
            binmean_PF_temp(isnan(binmean_PF_temp))=0;
            
            mean_PF=mean(binmean_PF_temp,2);
            mean_PF_withnoise=mean(binmean_temp,2);
            
            PF_center(PF_count)=find(mean_PF==max(mean_PF));
            PF_mean(PF_count,:)=mean_PF;
            PF_count=PF_count+1;
            placecell_id{f}=[placecell_id{f} i];
        end
    
    end
                   
end
%MEAN=[MEAN; mean_trans(placecell_id{f},:)];
end

[PF_sorted,sort_id]=sort(PF_center);
PF_mean_sorted=PF_mean(sort_id,:);
familiar_sort=sort_id;
f=figure;imagesc(PF_mean_sorted,[0,1.8]);
title('familiar')

%%


[behavior_filepaths, temp]=uigetfile('*.mat', 'Chose n files to load:','MultiSelect','on');
behavior_filepaths=sort(behavior_filepaths);
PF_center=[];
PF_mean=[];
PF_count=1;
for f=1:size(behavior_filepaths,2)
load([temp behavior_filepaths{f}]);
for i=1:size(sig_PFs,2)
    for j=1:size(sig_PFs,1)
        if ~isempty(sig_PFs{j,i})
            binmean_temp=sig_PFs_with_noise{j,i};
            binmean_PF_temp=sig_PFs{j,i};
            
            mean_PF=mean(binmean_PF_temp,2);
            mean_PF_withnoise=mean(binmean_temp,2);
            
            PF_center(PF_count)=find(mean_PF==max(mean_PF));
            PF_mean(PF_count,:)=mean_PF;
            PF_count=PF_count+1;

        end
    
    end
                   
end
end
[PF_sorted,sort_id]=sort(PF_center);
PF_mean_sorted=PF_mean(sort_id,:);

f=figure;imagesc(PF_mean_sorted,[0,1.8]);
title('novel')

MEAN=[];

for f=1:size(behavior_filepaths,2)
load([temp behavior_filepaths{f}]);
field_id=placecell_id{f};
MEAN=[MEAN; mean_trans(placecell_id{f},:)];

end
f=figure;imagesc(MEAN(familiar_sort,:),[0,1.8]);
title('f2n')