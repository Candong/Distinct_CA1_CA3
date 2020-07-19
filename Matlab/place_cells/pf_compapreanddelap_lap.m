clear all; close all;

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

% calculate the percentage of same pf in each plane
for i=1:3
    f_total(i)=length(f_first25_pfid{i})+length(f_late15_pfid{i})-length(f_potential_samepf{i});
    f_samepercent(i)=length(f_potential_samepf{i})/f_total(i);
      
end

figure;
labels={'common', '1st25only', 'late15only'};
for i=1:3
subplot(2,3,i)

f_1standlate=[rectl(length(f_potential_samepf{i})) rectl(length(f_first25_pfid{i})-length(f_potential_samepf{i})) rectl(length(f_late15_pfid{i})-length(f_potential_samepf{i}))];
pie(f_1standlate,labels)
title(['plane' num2str(i)]);

subplot(2,3,i+3)
bar([length(f_first25_pfid{i}) length(f_late15_pfid{i})])
xticks(1:2);
xticklabels({'1st','late'});

end



%% novel 1st 25 and late 15


for i=1:size(PF_info,2)
  if contains(PF_info(i).filepaths,'_n_')
      
      if contains(PF_info(i).filepaths,'1st25lap')
          
          if contains(PF_info(i).filepaths,'plain1')
          n_first25_pfid{1}=PF_info(i).place_cell_id;
          n_first25_id{1}=i;
          elseif contains(PF_info(i).filepaths,'plain2')
          n_first25_pfid{2}=PF_info(i).place_cell_id;
          n_first25_id{2}=i;          
          elseif contains(PF_info(i).filepaths,'plain3')
          n_first25_pfid{3}=PF_info(i).place_cell_id;
          n_first25_id{3}=i; 
          end
          
      elseif contains(PF_info(i).filepaths,'late15lap')
          
          if contains(PF_info(i).filepaths,'plain1')
          n_late15_pfid{1}=PF_info(i).place_cell_id;
          n_late15_id{1}=i;
          elseif contains(PF_info(i).filepaths,'plain2')
          n_late15_pfid{2}=PF_info(i).place_cell_id;
          n_late15_id{2}=i;
          elseif contains(PF_info(i).filepaths,'plain3')
          n_late15_pfid{3}=PF_info(i).place_cell_id;
          n_late15_id{3}=i;      
          end
          
      end
      
  end
    
end

for i=1:3
n_potential_samepf{i}=intersect(n_first25_pfid{i},n_late15_pfid{i});

end
% make sure all potential same pf is the same
for p=1:3
for i=n_potential_samepf{p} 
    if isempty(i)
      break
    end
load([temp behavior_filepaths{n_first25_id{p}}]);
    for j=1:size(sig_PFs,1)
        
        if ~isempty(sig_PFs{j,i})
       n_first25_pf(1)=PF_info(n_first25_id{p}).PF_start_bin(j,i);
       n_first25_pf(2)=PF_info(n_first25_id{p}).PF_end_bin(j,i); 
        end

    end
    
load([temp behavior_filepaths{n_late15_id{p}}]);    
    for j=1:size(sig_PFs,1)
        if ~isempty(sig_PFs{j,i}) 
       n_late15_pf(1)=PF_info(n_late15_id{p}).PF_start_bin(j,i);
       n_late15_pf(2)=PF_info(n_late15_id{p}).PF_end_bin(j,i);
        end

    end
pf_diff=sum(abs(n_first25_pf-n_late15_pf));

if pf_diff>4
    n_potential_samepf{p}(i)=[];
end
    
end
PF_info(n_first25_id{p}).same_pfid=n_potential_samepf{p};
end


% calculate the percentage of same pf in each plane
for i=1:3
    n_total(i)=length(n_first25_pfid{i})+length(n_late15_pfid{i})-length(n_potential_samepf{i});
    n_samepercent(i)=length(n_potential_samepf{i})/n_total(i);
      
end

figure;
labels={'common', '1st25only', 'late15only'};
for i=1:3
subplot(2,3,i)

n_1standlate=[rectl(length(n_potential_samepf{i})) rectl(length(n_first25_pfid{i})-length(n_potential_samepf{i})) rectl(length(n_late15_pfid{i})-length(n_potential_samepf{i}))];
pie(n_1standlate,labels)
title(['plane' num2str(i)]);

subplot(2,3,i+3)
bar([length(n_first25_pfid{i}) length(n_late15_pfid{i})])
xticks(1:2);
xticklabels({'1st','late'});

end


%% novel 1st 25 and 2nd 25

for i=1:size(PF_info,2)
  if contains(PF_info(i).filepaths,'_n_')
      
      if contains(PF_info(i).filepaths,'1st25lap')
          
          if contains(PF_info(i).filepaths,'plain1')
          n_first25_pfid{1}=PF_info(i).place_cell_id;
          n_first25_id{1}=i;
          elseif contains(PF_info(i).filepaths,'plain2')
          n_first25_pfid{2}=PF_info(i).place_cell_id;
          n_first25_id{2}=i;          
          elseif contains(PF_info(i).filepaths,'plain3')
          n_first25_pfid{3}=PF_info(i).place_cell_id;
          n_first25_id{3}=i; 
          end
          
      elseif contains(PF_info(i).filepaths,'2nd25lap')
          
          if contains(PF_info(i).filepaths,'plain1')
          n_2nd25_pfid{1}=PF_info(i).place_cell_id;
          n_2nd25_id{1}=i;
          elseif contains(PF_info(i).filepaths,'plain2')
          n_2nd25_pfid{2}=PF_info(i).place_cell_id;
          n_2nd25_id{2}=i;
          elseif contains(PF_info(i).filepaths,'plain3')
          n_2nd25_pfid{3}=PF_info(i).place_cell_id;
          n_2nd25_id{3}=i;      
          end
          
      end
      
  end
    
end

for i=1:3
n2_potential_samepf{i}=intersect(n_first25_pfid{i},n_2nd25_pfid{i});

end
% make sure all potential same pf is the same
for p=1:3
for i=n2_potential_samepf{p} 
    if isempty(i)
      break
    end
load([temp behavior_filepaths{n_first25_id{p}}]);
    for j=1:size(sig_PFs,1)
        
        if ~isempty(sig_PFs{j,i})
       n_first25_pf(1)=PF_info(n_first25_id{p}).PF_start_bin(j,i);
       n_first25_pf(2)=PF_info(n_first25_id{p}).PF_end_bin(j,i); 
        end

    end
    
load([temp behavior_filepaths{n_2nd25_id{p}}]);    
    for j=1:size(sig_PFs,1)
        if ~isempty(sig_PFs{j,i}) 
       n_2nd25_pf(1)=PF_info(n_2nd25_id{p}).PF_start_bin(j,i);
       n_2nd25_pf(2)=PF_info(n_2nd25_id{p}).PF_end_bin(j,i);
        end

    end
pf_diff=sum(abs(n_first25_pf-n_2nd25_pf));

if pf_diff>4
    n2_potential_samepf{p}(i)=[];
end
    
end
PF_info(n_2nd25_id{p}).same_pfid=n2_potential_samepf{p};
end


% calculate the percentage of same pf in each plane
for i=1:3
    n2_total(i)=length(n_first25_pfid{i})+length(n_2nd25_pfid{i})-length(n2_potential_samepf{i});
    n2_samepercent(i)=length(n2_potential_samepf{i})/n_total(i);
      
end
figure;
labels={'common', '1st25only', 'late15only'};
for i=1:3
subplot(2,3,i)

n2_1stand2nd=[rectl(length(n2_potential_samepf{i})) rectl(length(n_first25_pfid{i})-length(n2_potential_samepf{i})) rectl(length(n_2nd25_pfid{i})-length(n2_potential_samepf{i}))];
pie(n2_1stand2nd,labels)
title(['plane' num2str(i)]);

subplot(2,3,i+3)
bar([length(n_first25_pfid{i}) length(n_2nd25_pfid{i})])
xticks(1:2);
xticklabels({'1st','2nd'});

end

%% function for change negetive num to 0

function a=rectl(num)
if num<=0
    a=0;
else
    a=num;
end


end