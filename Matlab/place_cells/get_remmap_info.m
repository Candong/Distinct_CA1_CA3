clear all; close all;

[PFinfo_filepaths, temp]=uigetfile('*.mat', 'Chose info files to load:','MultiSelect','on');
load([temp PFinfo_filepaths]);
remap_info=struct;

for i=1:size(PF_info,2)
  if contains(PF_info(i).filepaths,'_f_')
      
     if contains(PF_info(i).filepaths,'late15lap')
          
          if contains(PF_info(i).filepaths,'plain1')
          remap_info(1).name='f plane1';
          remap_info(1).place_cellid=PF_info(i).combinePC;


          
          elseif contains(PF_info(i).filepaths,'plain2')
          remap_info(2).name='f plane2';
          remap_info(2).place_cellid=PF_info(i).combinePC;

          elseif contains(PF_info(i).filepaths,'plain3')
          remap_info(3).name='f plane3';
          remap_info(3).place_cellid=PF_info(i).combinePC;
          load(PF_info(i).filepaths);
          remap_info(3).meantrans=mean_trans;    
          end
     elseif  contains(PF_info(i).filepaths,'1st25lap')
          
          if contains(PF_info(i).filepaths,'plain1')
          remap_info(1).COM=PF_info(i).COM;
          remap_info(1).delay_lap=PF_info(i).delay_lap;
          load(PF_info(i).filepaths);
          remap_info(1).meantrans=mean_trans;
          elseif contains(PF_info(i).filepaths,'plain2')
          remap_info(2).COM=PF_info(i).COM;
          remap_info(2).delay_lap=PF_info(i).delay_lap;
          load(PF_info(i).filepaths);
          remap_info(2).meantrans=mean_trans;          
          elseif contains(PF_info(i).filepaths,'plain3')
          remap_info(3).COM=PF_info(i).COM; 
          remap_info(3).delay_lap=PF_info(i).delay_lap;
          load(PF_info(i).filepaths);
          remap_info(3).meantrans=mean_trans;
          
     end
  end
  end 
  if contains(PF_info(i).filepaths,'_n_')
      
     if contains(PF_info(i).filepaths,'late15lap')
          
          if contains(PF_info(i).filepaths,'plain1')
          remap_info(4).name='n plane1';
          remap_info(4).place_cellid=PF_info(i).combinePC;

          elseif contains(PF_info(i).filepaths,'plain2')
          remap_info(5).name='n plane2';
          remap_info(5).place_cellid=PF_info(i).combinePC;

          elseif contains(PF_info(i).filepaths,'plain3')
          remap_info(6).name='n plane3';
          remap_info(6).place_cellid=PF_info(i).combinePC;
    
          end
     elseif  contains(PF_info(i).filepaths,'1st25lap')
          
          if contains(PF_info(i).filepaths,'plain1')
          remap_info(4).COM=PF_info(i).COM;
          remap_info(4).delay_lap=PF_info(i).delay_lap;
          load(PF_info(i).filepaths);
          remap_info(4).meantrans=mean_trans;         
          elseif contains(PF_info(i).filepaths,'plain2')
          remap_info(5).COM=PF_info(i).COM;
          remap_info(5).delay_lap=PF_info(i).delay_lap;
          load(PF_info(i).filepaths);
          remap_info(5).meantrans=mean_trans;          
          elseif contains(PF_info(i).filepaths,'plain3')
          remap_info(6).COM=PF_info(i).COM;
          remap_info(6).delay_lap=PF_info(i).delay_lap;
          load(PF_info(i).filepaths);
          remap_info(6).meantrans=mean_trans;           
      end
     end

end

end
save([PF_info(1).filepaths(1:17) 'remap_info'],'remap_info');






