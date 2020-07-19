clear all; %close all;

[behavior_filepaths, temp]=uigetfile('*.mat', 'Chose f files to load:','MultiSelect','on');
PF_info=struct;
field_name={'plane1','plane2','plane3'};

PF_mean=[];
PF_count=1;
window=6;
threshold=3;
PF_start_lap=[];
range=1:30;
for f=1:size(behavior_filepaths,2)
load([temp behavior_filepaths{f}]);
place_cell_id=[];
for i=1:size(sig_PFs,2)  
    for j=1:size(sig_PFs,1)
        if ~isempty(sig_PFs{j,i})
            place_cell_id=[place_cell_id i];
        end

        end
    
end


PF_info(f).filepaths=behavior_filepaths{f};
PF_info(f).place_cell_id=place_cell_id; 
PF_info(f).PF_start_bin=PF_start_bins;
PF_info(f).PF_end_bin=PF_end_bins;

end


%% getnumPF matrix
name={'f1st25','n1st25','flate15','nlate15','n2nd25'};
num_PF=zeros(3,5);
for i=1:size(PF_info,2)
   if  contains(PF_info(i).filepaths,'_f_')
       
       if contains(PF_info(i).filepaths,'1st25lap')
           
           if contains(PF_info(i).filepaths,'plain1')
               num_PF(1,1)=length(PF_info(i).place_cell_id);
           elseif contains(PF_info(i).filepaths,'plain2')
               num_PF(2,1)=length(PF_info(i).place_cell_id);
           elseif contains(PF_info(i).filepaths,'plain3')
               num_PF(3,1)=length(PF_info(i).place_cell_id);
           end 
           
       elseif contains(PF_info(i).filepaths,'late15lap')
           
           if contains(PF_info(i).filepaths,'plain1')
               num_PF(1,3)=length(PF_info(i).place_cell_id);
           elseif contains(PF_info(i).filepaths,'plain2')
               num_PF(2,3)=length(PF_info(i).place_cell_id);
           elseif contains(PF_info(i).filepaths,'plain3')
               num_PF(3,3)=length(PF_info(i).place_cell_id);
           end  
           
       end
       
   elseif contains(PF_info(i).filepaths,'_n_')
       
       if contains(PF_info(i).filepaths,'1st25lap')
           
           if contains(PF_info(i).filepaths,'plain1')
               num_PF(1,2)=length(PF_info(i).place_cell_id);
           elseif contains(PF_info(i).filepaths,'plain2')
               num_PF(2,2)=length(PF_info(i).place_cell_id);
           elseif contains(PF_info(i).filepaths,'plain3')
               num_PF(3,2)=length(PF_info(i).place_cell_id);
           end 
           
       elseif contains(PF_info(i).filepaths,'late15lap')
           
           if contains(PF_info(i).filepaths,'plain1')
               num_PF(1,4)=length(PF_info(i).place_cell_id);
           elseif contains(PF_info(i).filepaths,'plain2')
               num_PF(2,4)=length(PF_info(i).place_cell_id);
           elseif contains(PF_info(i).filepaths,'plain3')
               num_PF(3,4)=length(PF_info(i).place_cell_id);
           end 
           
       elseif contains(PF_info(i).filepaths,'2nd25lap')
           
           if contains(PF_info(i).filepaths,'plain1')
               num_PF(1,5)=length(PF_info(i).place_cell_id);
           elseif contains(PF_info(i).filepaths,'plain2')
               num_PF(2,5)=length(PF_info(i).place_cell_id);
           elseif contains(PF_info(i).filepaths,'plain3')
               num_PF(3,5)=length(PF_info(i).place_cell_id);
           end  
           
       end
              
   end
   
        
end

%% compare different situation
%   f1st25 and late15 laps

for i=1:size(PF_info,2)
  if contains(PF_info(i).filepaths,'_f_')
      
      if contains(PF_info(i).filepaths,'1st25lap')
          
          if contains(PF_info(i).filepaths,'plain1')
          f_first25_pfid{1}=PF_info(i).place_cell_id;
          f_first25_id{1}=i;
          elseif contains(PF_info(i).filepaths,'plain2')
          f_first25_pfid{2}=PF_info(i).place_cell_id;
          f_first25_id{2}=i;          
          elseif contains(PF_info(i).filepaths,'plain3')
          f_first25_pfid{3}=PF_info(i).place_cell_id;
          f_first25_id{3}=i; 
          end
          
      else
          
          if contains(PF_info(i).filepaths,'plain1')
          f_late15_pfid{1}=PF_info(i).place_cell_id;
          f_late15_id{1}=i;
          elseif contains(PF_info(i).filepaths,'plain2')
          f_late15_pfid{2}=PF_info(i).place_cell_id;
          f_late15_id{2}=i;
          elseif contains(PF_info(i).filepaths,'plain3')
          f_late15_pfid{3}=PF_info(i).place_cell_id;
          f_late15_id{3}=i;      
          end
          
      end
      
  end
    
end

for i=1:3
f_potential_samepf{i}=intersect(f_first25_pfid{i},f_late15_pfid{i});

end
% make sure all potential same pf is the same
for p=1:3
for i=f_potential_samepf{p} 
load([temp behavior_filepaths{f_first25_id{p}}]);
    for j=1:size(sig_PFs,1)
        
        if ~isempty(sig_PFs{j,i})
       f_first25_pf(1)=PF_info(f_first25_id{p}).PF_start_bin(j,i);
       f_first25_pf(2)=PF_info(f_first25_id{p}).PF_end_bin(j,i); 
        end

    end
    
load([temp behavior_filepaths{f_late15_id{p}}]);    
    for j=1:size(sig_PFs,1)
        if ~isempty(sig_PFs{j,i}) 
       f_late15_pf(1)=PF_info(f_late15_id{p}).PF_start_bin(j,i);
       f_late15_pf(2)=PF_info(f_late15_id{p}).PF_end_bin(j,i);
        end

    end
pf_diff=sum(abs(f_first25_pf-f_late15_pf));

if pf_diff>4
    f_potential_samepf{p}(i)=[];
end
    
end
PF_info(f_first25_id{p}).same_pfid=f_potential_samepf{p};
end

%% find delay laps
window=6;
threshold=2;
range=1:25;
figure;
for plane=1:3
 behaviorfile=dir(['*ds_plain' num2str(plane) '.mat']);
 load(behaviorfile.name);
 cellsortfile=dir(['*Plain_' num2str(plane) '_*MotCor*' '.mat']);
 load(cellsortfile.name);
 neuron_id=combine2array(f_first25_pfid{plane},f_potential_samepf{plane}, f_late15_pfid{plane});

 %find start and end_bin for each place cell
 start_bin=[];end_bin=[];
 first25_id=setxor(f_first25_pfid{plane},f_potential_samepf{plane});
 late15_id=f_late15_pfid{plane};
 if ~isempty(first25_id)
 for k=first25_id
     load([temp behavior_filepaths{f_first25_id{plane}}]);
    for j=1:size(sig_PFs,1)
        
        if ~isempty(sig_PFs{j,k})
       start_bin=[start_bin PF_info(f_first25_id{plane}).PF_start_bin(j,k)];
       end_bin=[end_bin PF_info(f_first25_id{plane}).PF_end_bin(j,k)]; 
        end

    end   
 end
 end
 
  for k=late15_id
     load([temp behavior_filepaths{f_late15_id{plane}}]);
    for j=1:size(sig_PFs,1)
        
        if ~isempty(sig_PFs{j,k})
       start_bin=[start_bin PF_info(f_late15_id{plane}).PF_start_bin(j,k)];
       end_bin=[end_bin PF_info(f_late15_id{plane}).PF_end_bin(j,k)]; 
        end

    end   
  end
 
 PF_info(f_late15_id{plane}).combinePC=neuron_id;
 PF_info(f_late15_id{plane}).combinePC_start=start_bin;
 PF_info(f_late15_id{plane}).combinePC_end=end_bin;
 
 %analyze specific range:
 num_of_lap=1:25;
% startF=behavior.startframe(1);
% endF=length(behavior.ybinned);
startF=1;
endF=7000;
Fc3_DF=data.Fc3(startF:endF,:);
Fc2=data.Fc(startF:endF,:);
ybinned=behavior.ybinned(startF:endF);
% 
% empty_id=13610:14290-startF;
% Fc3_DF(empty_id,:)=[];
% Fc2(empty_id,:)=[];
% ybinned(empty_id,:)=[];



%%




ybinned_GoodBehav=ybinned;
Fc3_DF_GoodBehav=Fc3_DF;



%%%%%%%%%%%%%%%%%%%parameters%%%%%%%%%%%%%%%%%

ybinned=ybinned';
E=bwlabel(double_thresh(ybinned,0.0145,0.014)); %labels each traversal

%numneurons=size(Fc3_DF,2);
trackstart=min(ybinned_GoodBehav)+0.005; %track start location in quake units (+10 accounts for any noise in the track start location after teleportation)
trackend=max(ybinned_GoodBehav)-0.005; %track end location in quake units

numbins=50;

binsize=(trackend-trackstart)/numbins;%in Quake unit


%% correct E
    wrong_lap=0;
    for i=1:max(E)

            
        onpoint=find(E==i,1);
        offpoint=find(E==i,1,'last');    
%         end
        count=1;
        while max(ybinned(onpoint:offpoint))<0.638 
           count=count+1;
           if i~=max(E) 
             offpoint=find(E==i+1,1,'last');
             E(onpoint:offpoint)=i;
             E(offpoint+1:end)=E(offpoint+1:end)-1;
           else
               E(onpoint:offpoint)=0;
           end

        end
        
    end

    
    for i=1:max(E)

        if i==1
        onpoint=find(E==i,1);
        if ybinned(onpoint)>0.02
            E(E~=0)=E(E~=0)-1;
        end
        
        end

    end
E(E<min(num_of_lap)|E>max(num_of_lap))=0;
E=E-min(num_of_lap)+1;
mean_trans=[];
%%
count=1;
for ii = neuron_id
    
    %%%%%%%%now cut out transients from each lap and bin them%%%%%%%%%%%%%%%%%
    
    %seperate transients into laps
    cut_mat=NaN(2000,1000);
    cut_mat_ybinned=NaN(2000,1000);
    wrong_lap=0;
    for i=1:max(E)

        onpoint=find(E==i,1);
        offpoint=find(E==i,1,'last');    

        if max(ybinned(onpoint:offpoint))>0.6
        %cut out the transient
        cut_mat(1:(offpoint-onpoint)+1,i)=Fc3_DF(onpoint:offpoint,ii); %cut out activity from each lap
        %cut out the ybinned associated with the transient
        cut_mat_ybinned(1:(offpoint-onpoint)+1,i)=ybinned(onpoint:offpoint,1);
        %figure; hold on;plot(cut_mat_ybinned);
        wrong_lap=0;
        else
        wrong_lap=1;    
        end
        
    end
    for i=1:max(E)
        [cut_mat_ybinned_sorted(:,i), idx] = sort(cut_mat_ybinned(:,i),1);
        cut_mat_sorted(:,i)=cut_mat(idx,i);
    end
    
    %bin activity
    cut_mat_ybinned_sorted(cut_mat_ybinned_sorted==0)=nan;
    topEdge = trackend; % define limits
    botEdge = trackstart; % define limits
    binMean=[];
    for r=1:max(E)
        x = cut_mat_ybinned_sorted(:,r); %split into x and y
        y = cut_mat_sorted(:,r);
        binEdges = linspace(botEdge, topEdge, numbins+1);
        [h,whichBin] = histc(x, binEdges);
        for i = 1:numbins
            flagBinMembers = (whichBin == i);
            binMembers     = y(flagBinMembers);
            binMean(i,r)     = mean(binMembers);
        end
    end
    
%     mean_trans(ii,:)=mean(binMean,2);
%     figure;
%     imagesc(binMean');

    delay_lapid{plane}(count)=find_delaylap(binMean,start_bin(count),end_bin(count),window,threshold);
    
count=count+1;
    %close all
end %loops through each neuron
 
    PF_info(f_late15_id{plane}).delay_lap=delay_lapid{plane};
 
f_distribution=histc(delay_lapid{plane},range);
subplot(1,3,plane)
bar(f_distribution)
title(['familiar plane' num2str(plane)]);
ylim([0 25]);
 
 
end



%% function for change negetive num to 0

function a=rectl(num)
if num<=0
    a=0;
else
    a=num;
end


end


function newarray=combine2array(A,B,D)

C=setxor(A,B);
if size(C,1)>1
    C=C'
end
newarray=[C D];

end

