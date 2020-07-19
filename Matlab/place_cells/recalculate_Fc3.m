clear all;

[activity_filepaths, temp]=uigetfile('*.mat', 'Chose files for recalculate Fc3:','MultiSelect','on');

if ~isa(activity_filepaths,'cell') 
  activity_filepaths={activity_filepaths};  
end 

norm_window=500; % ms window used for f2fc normoliazation
end_frame=13000;

for f=1:size(activity_filepaths,2)
    
  load([temp activity_filepaths{f}]);
  
  Fc_1=FtoFc(data.F(1:end_frame,:),norm_window);
  Fc_2=FtoFc(data.F(end_frame+1:end,:),norm_window);
  
  Fc3_1=FctoFc3(Fc_1,upperbase(Fc_1,0));
  Fc3_2=FctoFc3(Fc_2,upperbase(Fc_2,0));
  con_Fc3=[Fc3_1; Fc3_2];
  
  data.Fc3_new=con_Fc3;
  
  if contains(activity_filepaths{f},'Plain_1')
    save('test1_new2','data');  
    
  elseif contains(activity_filepaths{f},'Plain_2')
    save('test2_new2','data'); 
    
  elseif contains(activity_filepaths{f},'Plain_3')
    save('test3_new2','data'); 
    
  end
  
  
end
