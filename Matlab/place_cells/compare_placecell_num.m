%clear all; close all;

[behavior_filepaths, temp]=uigetfile('*.mat', 'Chose files to load:','MultiSelect','on');
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


%% compare 1st25lap of f and novel
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

%% plot info
figure
bar(num_PF','stacked');
xticks(1:5);
xticklabels(name);
legend('plane1','plane2','plane3')



%% compare place field plot_each_heatmap.mplot_each_heatmap.m


