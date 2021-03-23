clear all; %close all;
% for now the slope is calculated by onstlap and later com differences, it
% could also include the earlier activities before 
% 
[pf_files, temp]=uigetfile('*.mat', 'Chose PF files to load:','MultiSelect','on');

f_count=1;
n_count=1;
for f=1:size(pf_files,2)
    if contains(pf_files{f},'_f_')
       f_pf_files{f_count}= pf_files{f};
       f_count=f_count+1;
    elseif contains(pf_files{f},'_n_')
        n_pf_files{n_count}= pf_files{f};
        n_count=n_count+1;
    
    end
end
f_pf_files=sort(f_pf_files);
n_pf_files=sort(n_pf_files);
%%
stepsize_com=5;
stepsize_skew=1;
onsetlap=1;
corr_thresh=0.5;
delta_skewness=[];
delta_COM=[];
slope=[];
sig=[];
skewness_slope=[];
skewness_sig=[];
r_square=[];
c_r_square=[];
count=0;
for f=1:size(f_pf_files,2)
load(f_pf_files{f});
f_sigPFs=sig_PFs;
f_meantrans=mean_trans;
load(n_pf_files{f});
n_sigPFs=sig_PFs;
n_meantrans=mean_trans;

[f_pf_num f_pf_id]=find(cellfun(@isempty,f_sigPFs)==0);
[n_pf_num n_pf_id]=find(cellfun(@isempty,n_sigPFs)==0);
common_PC_id=intersect(f_pf_id,n_pf_id);

cur_corr=diag(corr(f_meantrans(common_PC_id,:)', n_meantrans(common_PC_id,:)'));
highcor_id=common_PC_id(find(cur_corr>=corr_thresh));

for i=highcor_id'

    for j=1:size(f_sigPFs,1)
        if ~isempty(f_sigPFs{j,i})
            f_binmean{j}=f_sigPFs{j,i};
        end
        if ~isempty(n_sigPFs{j,i})
            n_binmean{j}=n_sigPFs{j,i};
        end        
    end
    
    if size(f_binmean,2)>1 | size(n_binmean,2)>1
        for m=1:size(f_binmean,2)
            for n=1:size(n_binmean,2)
                if (~isclipped(f_binmean{m})) & (~isclipped(n_binmean{n}))
                cur_pfcorr=corr(mean(f_binmean{m},2),mean(n_binmean{n},2));
                if cur_pfcorr>=corr_thresh
                    [f_skewness f_COM f_slope f_significance f_skew_slope skew_sig c_r s_r] = calculate_by_lap_skewness_COM(f_binmean{1},onsetlap,stepsize_com,1); % calculate variables for COM
                    [f_skewness nouseCOM nouse_slope nousesignificance f_skew_slope skew_sig nouser s_r] = calculate_by_lap_skewness_COM(f_binmean{1},onsetlap,stepsize_skew,1);% caluculate variables for skewness
                    [n_skewness n_COM] = calculate_by_lap_skewness_COM(n_binmean{1},onsetlap,stepsize_com,0);
                    [n_skewness nouseCOM] = calculate_by_lap_skewness_COM(n_binmean{1},onsetlap,stepsize_skew,0);
                    f_COM_removenan=f_COM(~isnan(f_COM));
                    n_COM_removenan=n_COM(~isnan(n_COM));
                    f_skewness_removenan=f_skewness(~isnan(f_skewness));
                    n_skewness_removenan=n_skewness(~isnan(n_skewness));
                    delta_skewness=[delta_skewness; f_skewness_removenan(1) f_skewness_removenan(end) n_skewness_removenan(1)];
                    delta_COM = [delta_COM;f_COM_removenan(1) f_COM_removenan(end) n_COM_removenan(1)];
                    slope=[slope f_slope];
                    sig=[sig f_significance];
                    skewness_slope=[skewness_slope f_skew_slope];
                    skewness_sig=[skewness_sig skew_sig];
                    r_square=[r_square s_r];
                    c_r_square=[c_r_square c_r];
                    count=count+1;
%                     figure; 
                    
%                     subplot(2,2,1)
%                     hold on
%                     plot(f_skewness,'o')
%                     plot(n_skewness,'o')
%                     subplot(2,2,2)
%                     hold on
%                     plot(f_COM)
%                     plot(n_COM)
%                     subplot(2,2,3)
%                     imagesc(f_binmean{m}')
%                     subplot(2,2,4)
%                     imagesc(n_binmean{n}')
%                     
%                     title(['corr=' num2str(cur_pfcorr)])
                end
                end
                
            end
        end
        
    else
        if (~isclipped(f_binmean{1})) & (~isclipped(n_binmean{1}))
            cur_pfcorr=corr(mean(f_binmean{1},2),mean(n_binmean{1},2));
            [f_skewness f_COM f_slope f_significance f_skew_slope skew_sig c_r s_r] = calculate_by_lap_skewness_COM(f_binmean{1},onsetlap,stepsize_com,1); % calculate variables for COM
            [f_skewness nouseCOM nouse_slope nousesignificance f_skew_slope skew_sig nouser s_r] = calculate_by_lap_skewness_COM(f_binmean{1},onsetlap,stepsize_skew,1);% caluculate variables for skewness
            [n_skewness n_COM] = calculate_by_lap_skewness_COM(n_binmean{1},onsetlap,stepsize_com,0);
            [n_skewness nouseCOM] = calculate_by_lap_skewness_COM(n_binmean{1},onsetlap,stepsize_skew,0);
            f_COM_removenan=f_COM(~isnan(f_COM));
            n_COM_removenan=n_COM(~isnan(n_COM));
            f_skewness_removenan=f_skewness(~isnan(f_skewness));
            n_skewness_removenan=n_skewness(~isnan(n_skewness));
            delta_skewness=[delta_skewness; f_skewness_removenan(1) f_skewness_removenan(end) n_skewness_removenan(1)];
            delta_COM = [delta_COM;f_COM_removenan(1) f_COM_removenan(end) n_COM_removenan(1)];
            slope=[slope f_slope];
            sig=[sig f_significance];
            skewness_slope=[skewness_slope f_skew_slope];
            skewness_sig=[skewness_sig skew_sig];
            r_square=[r_square s_r];
            c_r_square=[c_r_square c_r];
            count=count+1;
%             figure; 
%             subplot(1,2,1)
%             hold on
%             plot(f_skewness,'o')
%             plot(n_skewness,'o')
%             subplot(1,2,2)
%             hold on
%             plot(f_COM)
%             plot(n_COM)
%             
%             title(['corr=' num2str(cur_pfcorr)])
%         
        end
   
    end
    
    
end

% [slope all_start_lap COM_start COM_end COM_alllaps onset_deltaCOM all_deltaCOM pf_id]=caculate_shift_parameters(sig_PFs,window,onsetlap,step,fit_lapnum);
% f_slope=[f_slope slope];
% f_all_start_lap=[f_all_start_lap all_start_lap];
% f_COM_start=[f_COM_start COM_start];
% f_COM_end=[f_COM_end COM_end];
% f_COM_alllaps=[f_COM_alllaps COM_alllaps];
% f_pf_id{f}=pf_id;
% f_mean_trans{f}=mean_trans(pf_id,:);
% f_total_meantrans = [f_total_meantrans; mean_trans(pf_id,:)]; 
% f_all_deltaCOM {f}=all_deltaCOM;
% f_startlap{f}=all_start_lap;
% f_onset_deltaCOM=[f_onset_deltaCOM onset_deltaCOM];
end
%%
backward_id = slope<0 & sig;% & c_r_square>=0.05;
forward_id = slope >0 & sig;% & c_r_square>=0.05;
com_threshold=-0.1;
d0_com=delta_COM(:,1)-delta_COM(:,1);
d_com=delta_COM(:,2)-delta_COM(:,1);
d2_com=delta_COM(:,3)-delta_COM(:,1);

dd_com=delta_COM(:,3)-delta_COM(:,2);

day1_last=delta_COM(:,2);
day1_first=delta_COM(:,1);
day2_first=delta_COM(:,3);

%backward_id=backward_id' & d_com<0;
% d2_com(d_com>=com_threshold)=[];
% d0_com(d_com>=com_threshold)=[];
% dd_com(d_com>=com_threshold)=[];
% d_com(d_com>=com_threshold)=[];


figure;
plot([d0_com(backward_id)'; d_com(backward_id)'; d2_com(backward_id)'])
title('delta com')

figure;
subplot(1,3,1)
hold on
%histogram(delta_COM(:,3)-delta_COM(:,2),-49.5:1:49.5)
histogram(dd_com(backward_id),-13:1:13,'Normalization','Probability')
legend({'backward pf'})
title('delta com day2 fiirst - day1 last')

[p h]=signrank(dd_com(backward_id))
text(-2,0.1,[ 'signrank p=' num2str(p)]);

subplot(1,3,2)
plot(day1_last(backward_id),day2_first(backward_id),'o')
xlabel('day1 last lap')
ylabel('day2 first lap')
hold on;
plot([0 50],[0 50])

subplot(1,3,3)
plot(day1_first(backward_id),day2_first(backward_id),'o')
xlabel('day1 first lap')
ylabel('day2 first lap')
hold on;
plot([0 50],[0 50])

set(gcf, 'Position',  [100, 100, 1800, 500])

% figure;
% histogram(dd_com(backward_id)./d_com(backward_id))%,-5:0.1:5)
% 
% figure;
% histogram(dd_com(forward_id)./d_com(forward_id))

% figure;
% figure;
% histogram(dd_com(sig)./d_com(backward_id),-5:0.1:5)
% 
% figure;hold on;
% histogram(delta_COM(:,2)-delta_COM(:,1),-49.5:1:49.5)
% histogram(delta_COM(:,3)-delta_COM(:,1),-49.5:1:49.5)
% title('delta com')
%%
skew_backward_id = skewness_slope<0 & skewness_sig;%<0.05;% & r_square>=0.05;
skew_forward_id = skewness_slope >0 & skewness_sig;%<0.05;% & r_square>=0.05;



d0_skewness=delta_skewness(:,1)-delta_skewness(:,1);
d_skewness=delta_skewness(:,2)-delta_skewness(:,1);
d2_skewness=delta_skewness(:,3)-delta_skewness(:,1);

dd_skewness=delta_skewness(:,3)-delta_skewness(:,2);

%skew_backward_id=skew_backward_id' & d_skewness<0;
%dif_neg_id=d_skewness<0;
% d2_skewness(d_skewness>=skewness_threshold)=[];
% d0_skewness(d_skewness>=skewness_threshold)=[];
% dd_skewness(d_skewness>=skewness_threshold)=[];
% d_skewness(d_skewness>=skewness_threshold)=[];

figure;
plot([d0_skewness(skew_backward_id)'; d_skewness(skew_backward_id)'; d2_skewness(skew_backward_id)'])
title('delta skewness')

% figure;
% hold on
% %histogram(delta_COM(:,3)-delta_COM(:,2),-49.5:1:49.5)
% histogram(dd_skewness(skew_backward_id),-3:0.4:3,'Normalization','Probability')
% legend({'backward pf'})
% title('delta skewness')
% 
% [p h]=signrank(dd_skewness(skew_backward_id))
% text(-2,0.1,[ 'signrank p=' num2str(p)]);
% xlabel(' delta skewness (day2 -day1)')

day1_last=delta_skewness(:,2);
day1_first=delta_skewness(:,1);
day2_first=delta_skewness(:,3);

figure;
subplot(1,3,1)
hold on
%histogram(delta_COM(:,3)-delta_COM(:,2),-49.5:1:49.5)
histogram(dd_skewness(skew_backward_id),-3:0.4:3,'Normalization','Probability')
legend({'backward pf'})
title('delta skewness day2 fiirst - day1 last')

[p h]=signrank(dd_skewness(skew_backward_id))
text(-2,0.1,[ 'signrank p=' num2str(p)]);

subplot(1,3,2)
plot(day1_last(skew_backward_id),day2_first(skew_backward_id),'o')
xlabel('day1 last lap')
ylabel('day2 first lap')
hold on;
plot([-3 3],[-3 3])

subplot(1,3,3)
plot(day1_first(skew_backward_id),day2_first(skew_backward_id),'o')
xlabel('day1 first lap')
ylabel('day2 first lap')
hold on;
plot([-3 3],[-3 3])

set(gcf, 'Position',  [100, 100, 1800, 500])
%%
%%
% dif_neg_id=d_skewness<-0;
% 
% 
% figure;
% plot([d0_skewness(dif_neg_id)'; d_skewness(dif_neg_id)'; d2_skewness(dif_neg_id)'])
% title('delta skewness')
% 
% figure;
% hold on
% %histogram(delta_COM(:,3)-delta_COM(:,2),-49.5:1:49.5)
% histogram(dd_skewness(dif_neg_id),-3:0.2:3)
% legend({'backward pf'})
% title('delta skewness')
% 
% [p h]=signrank(dd_skewness(dif_neg_id))
% figure;
% histogram(dd_skewness(backward_id)./d_skewness(backward_id),-5:0.1:5)
% 
% figure;
% histogram(dd_skewness(forward_id)./d_skewness(forward_id))

% figure;hold on;
% %plot(delta_COM'-delta_COM(:,1)')
% histogram(delta_skewness(:,2)-delta_skewness(:,1),-9.5:0.1:9.5)
% histogram(delta_skewness(:,3)-delta_skewness(:,1),-9.5:0.1:9.5)
% title('delta skewness')
%%
% skewness_threshold=-0.2;
% d0_skewness=delta_skewness(:,1)-delta_skewness(:,1);
% d_skewness=delta_skewness(:,2)-delta_skewness(:,1);
% d2_skewness=delta_skewness(:,3)-delta_skewness(:,1);
% 
% dd_skewness=delta_skewness(:,3)-delta_skewness(:,2);
% 
% d2_skewness(d_skewness>=skewness_threshold)=[];
% d0_skewness(d_skewness>=skewness_threshold)=[];
% dd_skewness(d_skewness>=skewness_threshold)=[];
% d_skewness(d_skewness>=skewness_threshold)=[];
% 
% figure;
% plot([d0_skewness'; d_skewness'; d2_skewness'])
% title('delta skewness')
% 
% figure;hold on;
% %plot(delta_COM'-delta_COM(:,1)')
% %histogram(delta_skewness(:,3)-delta_skewness(:,2),-9.5:0.1:9.5)
% histogram(dd_skewness,-10:0.2:10)
% signrank(dd_skewness)
% legend({'backward pf'})
% title('delta skewness')