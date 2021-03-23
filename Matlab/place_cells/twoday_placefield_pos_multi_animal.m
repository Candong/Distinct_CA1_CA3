clear all; %close all;

%%
[PC_filepaths, temp]=sort(uigetfile('*.mat', 'Chose f & n files to load:','MultiSelect','on'));
planenum = input(['What are the planes you would like to analysis  >> ']);
summary_cor=[];
inst_summary=[];

f_count=1;
n_count=1;
window=6;
threshold=3;
max_lap=25;
unpass_count=0;


for f=1:size(PC_filepaths,2)
 if contains(PC_filepaths{f},'_f_')
     
        
    for p=planenum
        cur_delay_PF=[];
        cur_inst_PF=[];
        f_startlap=[];
        f_first_mean=[];
        f_second_mean=[];
        day1_last_mean=[];
    if contains(PC_filepaths{f},['plain' num2str(p)]) 

        load(PC_filepaths{f});
        f_file{p,f_count}=PC_filepaths{f};
        f_sig_PFs{p,f_count}=sig_PFs;
        f_sig_PFs_with_noise{p,f_count}=sig_PFs_with_noise;
        f_meantrans_withnoise{p,f_count}=mean_trans;
        f_PF_start_bins{p,f_count}=PF_start_bins;
        f_PF_end_bins{p,f_count}=PF_end_bins;
        
        
        [pf_num pf_id]=find(cellfun(@isempty,sig_PFs)==0);  
        pf_id=unique(pf_id);
        for i=1:length(pf_id)
        start=PF_start_bins(1,pf_id(i));
        PFend=PF_end_bins(1,pf_id(i));
        binmean=sig_PFs{1,pf_id(i)};
        noisy_binmean=sig_PFs_with_noise{1,pf_id(i)};
        max_lap=size(binmean,2);
        lap=find_delaylap(binmean,start,PFend,window,threshold,max_lap) ;
        f_startlap=[f_startlap lap];
        if lap>1
            cur_delay_PF=[cur_delay_PF pf_id(i)];
        elseif lap==1
            cur_inst_PF=[cur_inst_PF pf_id(i)];
        elseif lap<1
            unpass_count=unpass_count+1;
        end
        %if size(binmean,2)>=30
        half_lap=floor(max_lap/2);
        f_first_mean=[f_first_mean; mean(binmean(:,(1:half_lap)),2)'];
        f_second_mean=[f_second_mean; mean(binmean(:,end-half_lap+1:end),2)']; 
        day1_last_mean=[day1_last_mean;mean(noisy_binmean(:,end-half_lap+1:end),2)'];
        %end
        end
%         f_first_mean(f_startlap>10,:)=[];
%         f_second_mean(f_startlap>10,:)=[];
        
        F_first_mean{p,f_count}=f_first_mean;
        F_second_mean{p,f_count}=f_second_mean;
        day1_last{p,f_count}=day1_last_mean;
        
        f_delay_id{p,f_count}=cur_delay_PF;
        f_inst_id{p,f_count}=cur_inst_PF;
        f_start{p,f_count}=f_startlap;
        
        %if ~isempty(f_first_mean)
        f_count=f_count+1;  
        %end
    end
    end
 elseif contains(PC_filepaths{f},'_n_')

        
    for p=planenum
        cur_delay_PF=[];
        cur_inst_PF=[];
        n_startlap=[];
        n_first_mean=[];
        n_second_mean=[];
        day2_first_mean=[];
    if contains(PC_filepaths{f},['plain' num2str(p)]) 
        load(PC_filepaths{f});
        
        n_sig_PFs{p,n_count}=sig_PFs;
        n_sig_PFs_with_noise{p,n_count}=sig_PFs_with_noise;
        n_meantrans_withnoise{p,n_count}=mean_trans;
        n_PF_start_bins{p,n_count}=PF_start_bins;
        n_PF_end_bins{p,n_count}=PF_end_bins;
        
        [pf_num pf_id]=find(cellfun(@isempty,sig_PFs)==0);
        pf_id=unique(pf_id);
        for i=1:length(pf_id)
        start=PF_start_bins(1,pf_id(i));
        PFend=PF_end_bins(1,pf_id(i));
        binmean=sig_PFs{1,pf_id(i)};
        noisy_binmean=sig_PFs_with_noise{1,pf_id(i)};
        max_lap=size(binmean,2);
        lap=find_delaylap(binmean,start,PFend,window,threshold,max_lap) ;
        n_startlap=[n_startlap lap];
        if lap>1
            cur_delay_PF=[cur_delay_PF pf_id(i)];
        elseif lap==1
            cur_inst_PF=[cur_inst_PF pf_id(i)];
        end
        %if size(binmean,2)>=30
        half_lap=floor(max_lap/2);
        n_first_mean=[n_first_mean; mean(binmean(:,(1:half_lap)),2)'];
        n_second_mean=[n_second_mean; mean(binmean(:,end-half_lap+1:end),2)']; 
        day2_first_mean=[day2_first_mean; mean(noisy_binmean(:,1:half_lap),2)']; 
        %end
        
        end
%         n_first_mean(n_startlap>10,:)=[];
%         n_second_mean(n_startlap>10,:)=[];
        
        N_first_mean{p,n_count}=n_first_mean;
        N_second_mean{p,n_count}=n_second_mean;
        day2_first{p,n_count}=day2_first_mean;
        
        n_delay_id{p,n_count}=cur_delay_PF;
        n_inst_id{p,n_count}=cur_inst_PF;
        n_start{p,n_count}=n_startlap;
        
        %if ~isempty(n_first_mean)
           n_count=n_count+1; 
        %end
           
    end  
  end  
    
 end

end
%% 
total_cor=[];
exclude_total_cor=[];
common_total_cor=[];
total_2day_cor=[];
fn_total_cor=[];
inst_total_cor=[];
n_delay_total_cor=[];
f_inst_total_cor=[];
f_delay_total_cor=[];
total_f_mean=[];
total_n_mean=[];
common_total_start=[];
uncommon_total_start=[];
f_cor=[];
n_cor=[];
for i=1:f_count-1
    
for p=planenum
  
%   for i=1:size(f_sig_PFs{p},2)
if ~isempty(f_start{p,i}) && ~isempty(n_start{p,i})
[f_pf_num{p,i} f_pf_id{p,i}]=find(cellfun(@isempty,f_sig_PFs{p,i})==0);
[n_pf_num{p,i} n_pf_id{p,i}]=find(cellfun(@isempty,n_sig_PFs{p,i})==0);
f_pf_id{p,i}=unique(f_pf_id{p,i});
n_pf_id{p,i}=unique(n_pf_id{p,i});


cur_f_mean=f_meantrans_withnoise{p,i};
cur_n_mean=n_meantrans_withnoise{p,i};
cur_d1_mean=day1_last{p,i};
cur_d2_mean=day2_first{p,i};

total_f_mean=[total_f_mean; cur_f_mean];
total_n_mean=[total_n_mean; cur_n_mean];

f_n_PCid{p,i}=union(f_pf_id{p,i},n_pf_id{p,i});

transient_cor{p,i}=diag(corr(cur_f_mean(f_n_PCid{p,i},:)',cur_n_mean(f_n_PCid{p,i},:)'));

if ~isempty(F_first_mean{p,i}) && ~isempty(F_second_mean{p,i})
cur_f_cor=diag(corr(F_first_mean{p,i}',F_second_mean{p,i}'));
end
if ~isempty(N_first_mean{p,i}) && ~isempty(N_second_mean{p,i})

cur_n_cor=diag(corr(N_first_mean{p,i}',N_second_mean{p,i}'));
end

      
total_cor=[total_cor transient_cor{p,i}'];   
f_cor=[f_cor cur_f_cor'];
n_cor=[n_cor cur_n_cor'];

exclude_PC_id{p,i}=setxor(f_pf_id{p,i},n_pf_id{p,i});
exclude_cor{p,i}=diag(corr(cur_f_mean(exclude_PC_id{p,i},:)',cur_n_mean(exclude_PC_id{p,i},:)'));
exclude_total_cor=[exclude_total_cor exclude_cor{p,i}'];    


common_PC_id{p,i}=intersect(f_pf_id{p,i},n_pf_id{p,i});
day1_include_id=ismember(f_pf_id{p,i},common_PC_id{p,i});
day2_include_id=ismember(n_pf_id{p,i},common_PC_id{p,i});

if ~isempty(common_PC_id{p,i})
common_cor{p,i}=diag(corr(cur_f_mean(common_PC_id{p,i},:)',cur_n_mean(common_PC_id{p,i},:)'));
day1_day2_cor{p,i}=diag(corr(cur_d1_mean(day1_include_id,:)',cur_d2_mean(day2_include_id,:)'));
common_start=n_start{p,i}(ismember(n_pf_id{p,i},common_PC_id{p,i}));
uncommon_start=n_start{p,i}(~ismember(n_pf_id{p,i},common_PC_id{p,i}));
common_total_cor=[common_total_cor common_cor{p,i}']; 
common_total_start=[common_total_start common_start];
uncommon_total_start=[uncommon_total_start uncommon_start];
total_2day_cor=[total_2day_cor day1_day2_cor{p,i}'];

end
 
% common_PC_id{p,i}=intersect(f_pf_id{p,i},n_pf_id{p,i});
% if ~isempty(common_PC_id{p,i})
% common_cor{p,i}=diag(corr(cur_f_mean(common_PC_id{p,i},:)',cur_n_mean(common_PC_id{p,i},:)'));
% common_total_cor=[common_total_cor common_cor{p,i}']; 
% end


if ~isempty(f_pf_id{p,i})
fn_cor{p,i}=diag(corr(cur_f_mean(f_pf_id{p,i},:)',cur_n_mean(f_pf_id{p,i},:)'));
fn_total_cor=[fn_total_cor fn_cor{p,i}']; 

end

if ~isempty(n_inst_id{p,i})
    inst_cor{p,i}=diag(corr(cur_f_mean(n_inst_id{p,i},:)',cur_n_mean(n_inst_id{p,i},:)'));
    inst_total_cor=[inst_total_cor inst_cor{p,i}'];
    
end

if ~isempty(n_delay_id{p,i})
    n_delay_cor{p,i}=diag(corr(cur_f_mean(n_delay_id{p,i},:)',cur_n_mean(n_delay_id{p,i},:)'));
    n_delay_total_cor=[n_delay_total_cor n_delay_cor{p,i}'];
    
end

if ~isempty(f_inst_id{p,i})
    f_inst_cor{p,i}=diag(corr(cur_n_mean(f_inst_id{p,i},:)',cur_f_mean(f_inst_id{p,i},:)'));
    f_inst_total_cor=[f_inst_total_cor f_inst_cor{p,i}'];
    
end

if ~isempty(f_delay_id{p,i})
    f_delay_cor{p,i}=diag(corr(cur_n_mean(f_delay_id{p,i},:)',cur_f_mean(f_delay_id{p,i},:)'));
    f_delay_total_cor=[f_delay_total_cor f_delay_cor{p,i}'];
    
end
    
end
end
%summary_cor=[summary_cor fn_total_cor];
%inst_summary=[inst_summary inst_total_cor];
end
%%
figure;
total_cor(isnan(total_cor))=[];
f_cor(isnan(f_cor))=[];
n_cor(isnan(n_cor))=[];
f_cor_sem= std(f_cor)/sqrt(length(f_cor));
total_cor_sem= std(total_cor)/sqrt(length(total_cor));


subplot(1,2,1)
bar([1 2],[mean(f_cor) mean(total_cor)]);
hold on;
errorbar([1 2],[mean(f_cor) mean(total_cor)],[f_cor_sem total_cor_sem],'.','LineWidth',2);
subplot(1,2,2)
hold on;
plot([1],f_cor,'o')
plot([2],total_cor,'o')
plot([3],n_cor,'o')

line([1.9 2.1],[mean(total_cor) mean(total_cor)],'Linewidth',2);
line([0.9 1.1],[mean(f_cor) mean(f_cor)],'Linewidth',2);
line([2.9 3.1],[median(n_cor) median(n_cor)],'Linewidth',2);
line([2.9 3.1],[mean(n_cor) mean(n_cor)],'Color','red','Linewidth',2);


g1 = repmat({'F'},length(f_cor),1);
g2 = repmat({'N'},length(n_cor),1);
g3 = repmat({'FN'},length(total_cor),1);
g=[g1;g2;g3];
figure
boxplot([f_cor,n_cor,total_cor],g)
ylim([-1,1])


range=-1:0.1:1;
  figure;
histogram(f_cor,range)
hold on;
histogram(n_cor,range)
legend({'f','n'})
title('f & n 1st last 10 lap correlation')
%save('ca3_cor_halflap','f_cor','total_cor','n_cor');
%plot([1],mean(total_cor),'or','LineWidth',3)
% figure;
% plot([1],exclude_total_cor,'o')
% hold on;plot([1],mean(exclude_total_cor),'*')
% 
% 
% figure;
% plot([1],common_total_cor,'o')
% hold on;plot([1],mean(common_total_cor),'*')
% 
% 
% figure;
% plot([1],fn_total_cor,'o')
% hold on;plot([1],mean(fn_total_cor),'*')
%%
range=-1:0.1:1;
[fn_dist edges]=histcounts(fn_total_cor,range);
figure;
subplot(1,2,1)
histogram(common_total_cor,range)
subplot(1,2,2)
histogram(common_total_cor,range,'Normalization','probability')
ylim([0,0.3])
figure;
subplot(1,2,1)
histogram(total_2day_cor,range)
subplot(1,2,2)
histogram(total_2day_cor,range,'Normalization','probability')
ylim([0,0.35])

figure;
histogram(common_total_start(common_total_cor>=0.5),1:30)%,'Normalization','probability')

hold on;
histogram(common_total_start(common_total_cor<0.5),1:30)%,'Normalization','probability')
legend({'high cor','low cor'})
xlabel('PF onset')

%save('ca1_2day_cor','total_2day_cor');
%%
%save('CA1_n2day_common_info','common_total_cor','common_total_start');
common_total_start(isnan(total_2day_cor))=[];
total_2day_cor(isnan(total_2day_cor))=[];
remove_late=find(common_total_start>=30);
total_2day_cor(remove_late)=[];
common_total_start(remove_late)=[];
%%
figure;
[hico_startlap_dist edges]=histcounts(common_total_start(total_2day_cor>=0.5),1:30);
[lowco_startlap_dist edges]=histcounts(common_total_start(total_2day_cor<0.5),1:30);
bar(hico_startlap_dist/length(total_2day_cor),'BarWidth', 1,'FaceAlpha',0.6)
hold on;
bar(lowco_startlap_dist/length(total_2day_cor),'BarWidth', 1,'FaceAlpha',0.6)
legend({'high cor','low cor'})
ylabel('percentage')
xlabel('PF onset')
%%
consistent_start=common_total_start(total_2day_cor>=0.5);
new_start=[uncommon_total_start common_total_start(total_2day_cor<0.5)];
consistent_start(consistent_start>30)=[];
new_start(new_start>30)=[];
figure;
histogram(consistent_start,1:30,'Normalization','Probability');
title('CA1 nday2 consistent PF onset lap')
figure;
histogram(new_start,1:30,'Normalization','Probability');
%legend('consistent','newly formed')
title('CA1 nday2  newly formed PF onset lap')

%%
figure
cdfplot(consistent_start);
hold on; cdfplot(new_start);
xlim([0 30])
[p h stats]=ranksum(consistent_start,new_start)
title(['CA3 onsetlap wilcoxon test,p= ' num2str(p)]);
xlabel('PC onset lap')
legend('consistent','newly formed')
grid off

%%
figure
cdfplot(common_total_start(total_2day_cor>=0.5));
hold on; cdfplot(common_total_start(total_2day_cor<0.5));
xlim([0 30])
p=ranksum(common_total_start(total_2day_cor>=0.5),common_total_start(total_2day_cor<0.5));
title(['CA1 wilcoxon test,p= ' num2str(p)]);
xlabel('PC onset lap')
legend({'high cor','low cor'})
grid off

% figure;bar(fn_dist)
% title('fn correlation')

% figure;
% histogram(fn_total
% ylim([0 65])

% figure;
% % histogram(f_inst_total_cor,range)
% n_delay_dist=histc(f_inst_total_cor,range);
% figure;bar(n_delay_dist)
% title('f inst in n correlation')
% ylim([0 65])

% figure;
% % histogram(f_delay_total_cor,range)
% title('f delay in n correlation')
% ylim([0 65])


% figure;
% histogram(inst_total_cor,range)
n_inst_dist=histc(inst_total_cor,range);
% figure;bar(n_inst_dist)
% title('n inst in f correlation')
% ylim([0 65])

% figure;
% histogram(n_delay_total_cor,range)

% ylim([0 65])

n_delay_dist=histc(n_delay_total_cor,range);
% figure;bar(n_delay_dist);
% title('n delay in f correlation')
%ylim([0 65])
%%  
for i=1:size(f_inst_cor,2)
    [cur_instdist edges]=histcounts(f_inst_cor{i},range);
    f_instdist(i,:)=cur_instdist;
    [cur_delaydist edges]=histcounts(f_delay_cor{i},range);
    f_delaydist(i,:)=cur_delaydist;    
    
    
end

inst_percent=sum(f_instdist(:,end-1:end),2)./sum(f_instdist,2);
delay_percent=sum(f_delaydist(:,end-1:end),2)./sum(f_delaydist,2);
% [meanCOM cur_SP]=meanCOMandSP(PF_denoise,start_b,end_b,numbins);

%%
remove_id=find(mean(total_f_mean,2)>5);
total_f_mean(remove_id,:)=[];
total_n_mean(remove_id,:)=[];

pop_cor_bybin=corr(total_f_mean,total_n_mean);
% figure;
% imagesc(flip(pop_cor_bybin));
% title('pop corelation by bin')
% hold on;plot(1:50,50:-1:1,'LineWidth',2)
% colormap jet
%%
function [meanCOM SP]=meanCOMandSP(binmean,start_b,end_b,numbins)
binM=1:numbins;
COM=sum(binmean(:,start_b:end_b).*binM(start_b:end_b),2)./sum(binmean(:,start_b:end_b),2);
COM(isnan(COM))=0;
A=max(binmean(:,start_b:end_b),[],2);
meanCOM=sum(A.*COM)/sum(A);

SP=1/(sqrt(sum(A.*((COM-meanCOM).^2))/sum(A)));
end 
 



