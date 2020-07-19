clear all
close all

%%
[behavior_filepaths, temp]=uigetfile('*.mat', 'Chose f and n files to load:','MultiSelect','on');
PF_info=struct;
field_name={'plane1','plane2','plane3'};

PF_mean=[];
PF_count=1;

PF_start_lap=[];

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