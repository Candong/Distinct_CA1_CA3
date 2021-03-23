clear all; close all;

%%
f1_v=struct;n1_v=struct;
f2_v=struct;n2_v=struct;
nday2_v=struct;
%% 
[behavior_path_day1, temp]=uigetfile('*.mat', 'Chose switch day1 behavior files to load:','MultiSelect','on');
%behavior_filepaths=sort(behavior);
%% seperate laps by reward
if behavior_path_day1~=0
load([temp behavior_path_day1]);

if behavior.startframe(1)<500
    f_start=behavior.startframe(1);


else
    f_start=1;
    
end

binnum=50;
track_length=300;
sample_rate=10.4;

y_min=min(behavior.ybinned);
y_max=max(behavior.ybinned);
edges=linspace(y_min,y_max,binnum+1);
y_belong=discretize(behavior.ybinned,edges);

f_end=floor(behavior.startframe(end)/1000)*1000;
n_start=behavior.startframe(end);
n_end=length(behavior.ybinned);

f_y=behavior.ybinned(f_start:f_end);
f_reward=behavior.reward(f_start:f_end);
f_velocity=behavior.velocity(f_start:f_end);
f_y_belong=y_belong(f_start:f_end);

f1_v=find_velocity(f_y,f_reward,f_velocity,f_y_belong,sample_rate,track_length,binnum);

n_y=behavior.ybinned(n_start:n_end);
n_reward=behavior.reward(n_start:n_end);
n_velocity=behavior.velocity(n_start:n_end);
n_y_belong=y_belong(n_start:n_end);
n1_v=find_velocity(n_y,n_reward,n_velocity,n_y_belong,sample_rate,track_length,binnum);

% figure;hold on;
% plot(f_v.bylap);
% plot(n_v.bylap)
end
%%

[behavior_path_day2, temp]=uigetfile('*.mat', 'Chose switch day2 behavior files to load:','MultiSelect','on');

%% seperate laps by reward
if behavior_path_day2~=0
load([temp behavior_path_day2]);

if behavior.startframe(1)<500
    f_start=behavior.startframe(1);
    


else
    f_start=1;
    
end

binnum=50;
track_length=300;
sample_rate=10.4;

y_min=min(behavior.ybinned);
y_max=max(behavior.ybinned);
edges=linspace(y_min,y_max,binnum+1);
y_belong=discretize(behavior.ybinned,edges);

f_end=floor(behavior.startframe(end)/1000)*1000;
n_start=behavior.startframe(end);
n_end=length(behavior.ybinned);

f_y=behavior.ybinned(f_start:f_end);
f_reward=behavior.reward(f_start:f_end);
f_velocity=behavior.velocity(f_start:f_end);
f_y_belong=y_belong(f_start:f_end);

f2_v=find_velocity(f_y,f_reward,f_velocity,f_y_belong,sample_rate,track_length,binnum);

n_y=behavior.ybinned(n_start:n_end);
n_reward=behavior.reward(n_start:n_end);
n_velocity=behavior.velocity(n_start:n_end);
n_y_belong=y_belong(n_start:n_end);
n2_v=find_velocity(n_y,n_reward,n_velocity,n_y_belong,sample_rate,track_length,binnum);

end






%%
[behavior_path_nday2, temp]=uigetfile('*.mat', 'Chose novel_day2 behavior files to load:','MultiSelect','on');

if behavior_path_nday2~=0
load([temp behavior_path_nday2]);
% 
% if behavior.startframe(1)<500
%     f_start=behavior.startframe(1);
% 
% 
% else
%     f_start=1;
%     
% end
% 
% binnum=50;
% track_length=300;
% sample_rate=10.4;

% y_min=min(behavior.ybinned);
% y_max=max(behavior.ybinned);
% edges=linspace(y_min,y_max,binnum+1);
y_belong=discretize(behavior.ybinned,edges);

% f_end=floor(behavior.startframe(end)/1000)*1000;
if behavior.startframe(end)>7000
    n2_start=behavior.startframe(end);
    n2_end=length(behavior.ybinned);
else 
    n2_start=input('What is the startfrme of nday2? ');
    n2_end=length(behavior.ybinned);
end

n2_y=behavior.ybinned(n2_start:n2_end);
n2_reward=behavior.reward(n2_start:n2_end);
n2_velocity=behavior.velocity(n2_start:n2_end);
n2_y_belong=y_belong(n2_start:n2_end);

nday2_v=find_velocity(n2_y,n2_reward,n2_velocity,n2_y_belong,sample_rate,track_length,binnum);

end
%%
figure;
subplot(2,1,1)
hold on;
if behavior_path_day1~=0
plot(f1_v.bylap_nofreeze);
plot(n1_v.bylap_nofreeze);
legend({'f1','n1'},'Location','bestoutside')
title('velocity by lap remove freezeing (cm/s) ')
end

subplot(2,1,2)
hold on;
if behavior_path_day2~=0
plot(f2_v.bylap_nofreeze);
plot(n2_v.bylap_nofreeze);
end

if behavior_path_nday2~=0
plot(nday2_v.bylap_nofreeze);
end
legend({'f2','n2','nday2'},'Location','bestoutside')
title('velocity by lap remove freezeing (cm/s) ')

figure;
subplot(2,1,1)
hold on;
if behavior_path_day1~=0
plot(f1_v.bylap);
plot(n1_v.bylap);
legend({'f1','n1'},'Location','bestoutside')
title('velocity by lap with freezeing (cm/s)')
end
subplot(2,1,2)
hold on;
if behavior_path_day2~=0
plot(f2_v.bylap);
plot(n2_v.bylap);
end
if behavior_path_nday2~=0
plot(nday2_v.bylap);
end
legend({'f2','n2','nday2'},'Location','bestoutside')
title('velocity by lap with freezeing (cm/s)')

% figure;
% subplot(2,1,1)
% hold on;
% if behavior_path_day1~=0
% plot(f1_v.freeze_by_lap./sample_rate);
% plot(n1_v.freeze_by_lap./sample_rate);
% legend({'f1','n1'},'Location','bestoutside')
% title('freeze by lap (sec)')
% end
% subplot(2,1,2)
% hold on;
% if behavior_path_day2~=0
% plot(f2_v.freeze_by_lap./sample_rate);
% plot(n2_v.freeze_by_lap./sample_rate);
% end
% if behavior_path_nday2~=0
% plot(nday2_v.freeze_by_lap./sample_rate);
% end
% legend({'f2','n2','nday2'},'Location','bestoutside')
% title('freeze by lap (sec)')

save([behavior_path_day1(15:end-4) 'velocity_info_new'],'f1_v','n1_v','f2_v','n2_v','nday2_v');
%%
function v_info=find_velocity(ybinned,reward,velocity,y_belong,sample_rate,track_length,binnum)

    reward_id=find(reward>8);
    reward_id(find(diff(reward_id)==1))=[];
    v_bylap=[];
    v_bylap_nofreeze=[];
    freeze_bylap=[];
    figure;hold on;
    y_binmatrix=ones(length(reward_id),binnum).*(track_length/binnum);
    y_bin_N=zeros(length(reward_id),binnum);
    for i=1:length(reward_id)
        if i==1
            
            cur_v=track_length/(length(1:reward_id(i))/sample_rate);
            %plot(ybinned(1:reward_id(i)));
            cur_freeze=find(velocity(1:reward_id(i))<=0.005);
            cur_y=ybinned(1:reward_id(i));
            cur_y(cur_freeze)=[];
            cur_v_nofreeze=track_length/(length(cur_y)/sample_rate);
            plot(cur_y);
            
            cur_y_belong=y_belong(1:reward_id(i));
            cur_y_belong(cur_freeze)=[];
            
            for j=1:binnum
            y_bin_N(i,j)=sum(cur_y_belong==j);   
                
                
            end
            
                
        else
            cur_v=track_length/((length(reward_id(i-1):reward_id(i))/sample_rate)-1.5);
            %plot(ybinned(reward_id(i-1):reward_id(i)));
            cur_freeze=find(velocity((reward_id(i-1)+ceil(1.5*sample_rate)):reward_id(i))<=0.005);
            cur_y=ybinned((reward_id(i-1)+ceil(1.5*sample_rate)):reward_id(i));
            cur_y(cur_freeze)=[];
            cur_v_nofreeze=track_length/(length(cur_y)/sample_rate);
            plot(cur_y); 
            
            cur_y_belong=y_belong(reward_id(i-1):reward_id(i));
            cur_y_belong(cur_freeze)=[];
            
            for j=1:binnum
            y_bin_N(i,j)=sum(cur_y_belong==j);   
                
            end            
            
        end
        v_bylap=[v_bylap cur_v];
        v_bylap_nofreeze=[v_bylap_nofreeze cur_v_nofreeze];
        freeze_bylap=[freeze_bylap length(cur_freeze)];
       
    end
    v_info.bylap=v_bylap;
    v_info.bylap_nofreeze=v_bylap_nofreeze;
    v_info.freeze_by_lap=freeze_bylap;
    v_info.v_bybinedlap=y_binmatrix./(y_bin_N/sample_rate);
    v_info.sample_rate=sample_rate;
    v_info.binnum=binnum;
    v_info.tracklength=track_length;

end

