clear all;
close all;

prompt = 'how many plains does this dataset contain? ';
plain_num = input(prompt);
[behavior_filepaths, temp]=uigetfile('*.mat', 'Chose behavior files to load:','MultiSelect','on');
for i=1:plain_num;
fname{i}=[behavior_filepaths(1:end-4) '_plain' num2str(i)]; %[behavior_filepaths(1:end-4) '_plain2'] };

end

load([temp behavior_filepaths]);


prompt = 'Is this a remapping data file? 1=yes 0=no? ';
x = input(prompt);

% session.beh_data.Y_pos=[session.beh_data.Y_pos ;0 ;0];
% session.beh_data.reward=[session.beh_data.reward ;0 ;0];
% session.beh_data.t=[session.beh_data.t 0 0];
% session.beh_data.velocity=[session.beh_data.velocity; 0 ;0];
% session.beh_data.lick=[session.beh_data.lick; 0 ;0];

y_pos=reshape(session.beh_data.Y_pos,plain_num,[]);
reward=reshape(session.beh_data.reward,plain_num,[]);
t=reshape(session.beh_data.t,plain_num,[]);
velocity=reshape(session.beh_data.velocity,plain_num,[]);
lick=reshape(session.beh_data.lick,plain_num,[]);

figure;hold on
for i=1:plain_num
behavior=struct;

behavior.ybinned=y_pos(i,:);

behavior.reward=max(reward);

behavior.t=t(i,:);

behavior.velocity=velocity(i,:);

behavior.lick=lick(i,:);

if x==0
    behavior.startframe=1;

    
else
    id=find((diff(behavior.ybinned)<-0.03 &diff(behavior.ybinned)>-0.4)==1)+1;
    behavior.startframe=id(find(behavior.ybinned(id-5)<0.5));
  
end

%plot(behavior.ybinned);
hold on
% 
% if i~=1
%  behavior.ybinned(end)=[];
%  behavior.reward(end)=[];
%  behavior.t(end)=[];
%  behavior.velocity(end)=[];    
% end
save(fname{i},'behavior')
end
figure;
plot(y_pos')
hold on;
plot(reward')


mean_behavior=struct;
mean_behavior.ybinned=y_pos;

mean_behavior.reward=max(reward);

mean_behavior.t=t;

mean_behavior.velocity=velocity;

mean_behavior.lick=lick;

if x==0
    mean_behavior.startframe=1;

    
else
    id=find((diff(mean_behavior.ybinned)<-0.03 &diff(mean_behavior.ybinned)>-0.4)==1)+1;
    mean_behavior.startframe=id(find(mean_behavior.ybinned(id-5)<0.5));
  
end



save([behavior_filepaths(1:end-4) '_mean'],'mean_behavior')