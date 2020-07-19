%calculate persantage
clear all
%close all
% 
% %load files
[behavior, temp]=uigetfile('*.mat', 'Chose f and n files to load:','MultiSelect','on');
file_count=1;
for f=1:size(behavior,2)
    if contains(behavior{f},'_f_')
       behavior_filepaths{file_count}= behavior{f};
       file_count=file_count+1;
    
    end
end

behavior_filepaths=sort(behavior_filepaths);

last_file=behavior_filepaths{1};
animal_count=1;
group_filepath{1}{1}=behavior_filepaths{1};
group_count=1;
for f=2:size(behavior_filepaths,2)
    %load(behavior_filepaths{f}); 
    cur_file=behavior_filepaths{f};
    if sum(cur_file(1:25)==last_file(1:25))==25
        group_count=group_count+1;
        group_filepath{animal_count}{group_count}=cur_file;

    else
        animal_count=animal_count+1;
        group_count=1;
        group_filepath{animal_count}{group_count}=cur_file;
        last_file=cur_file;
               
    end
    
end

f_neuron_num=[];
f_PC_num=[];
for i=1:size(group_filepath,2)
    neuron_num=0;
    PC_num=0;
    
    for j=1:size(group_filepath{i},2)
        load(group_filepath{i}{j})
        neuron_num=neuron_num+size(mean_trans,1);
        PC_num=PC_num+sum(~isnan(number_of_PFs));
        
    end
    f_neuron_num=[f_neuron_num neuron_num];
    f_PC_num=[f_PC_num PC_num];
    
end
f_percent=f_PC_num./f_neuron_num;

%%
[behavior, temp]=uigetfile('*.mat', 'Chose f and n files to load:','MultiSelect','on');
file_count=1;
for f=1:size(behavior,2)
    if contains(behavior{f},'_n_')
       behavior_filepaths{file_count}= behavior{f};
       file_count=file_count+1;
    
    end
end

behavior_filepaths=sort(behavior_filepaths);

last_file=behavior_filepaths{1};
animal_count=1;
group_filepath{1}{1}=behavior_filepaths{1};
group_count=1;
for f=2:size(behavior_filepaths,2)
    %load(behavior_filepaths{f}); 
    cur_file=behavior_filepaths{f};
    if sum(cur_file(1:25)==last_file(1:25))==25
        group_count=group_count+1;
        group_filepath{animal_count}{group_count}=cur_file;

    else
        animal_count=animal_count+1;
        group_count=1;
        group_filepath{animal_count}{group_count}=cur_file;
        last_file=cur_file;
               
    end
    
end

n_neuron_num=[];
n_PC_num=[];
for i=1:size(group_filepath,2)
    neuron_num=0;
    PC_num=0;
    
    for j=1:size(group_filepath{i},2)
        load(group_filepath{i}{j})
        neuron_num=neuron_num+size(mean_trans,1);
        PC_num=PC_num+sum(~isnan(number_of_PFs));
        
    end
    n_neuron_num=[n_neuron_num neuron_num];
    n_PC_num=[n_PC_num PC_num];
    
end
n_percent=n_PC_num./n_neuron_num;

%% 
figure;

bar([mean(f_percent) mean(n_percent)]/mean(f_percent))
hold on
% plot(1,f_percent,'o')
% plot(2,n_percent,'o')

plot([1 2],[f_percent' n_percent']./f_percent','o-')
%ylim([0 2])

figure;

bar([mean(f_percent); mean(n_percent)])
hold on
% plot(1,f_percent,'o')
% plot(2,n_percent,'o')

plot([1 2],[f_percent' n_percent'],'o-')