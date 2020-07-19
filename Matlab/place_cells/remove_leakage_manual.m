clear all;
close all;
%%
dist_threshold=45; % could be twiked with CA1 data
dist_step1=35;
dist_step2=25;
dist_step3=15;
base_threshold=0.92;
sim_threshold1=0.9;
sim_threshold2=0.92;
sim_threshold3=0.96;
sim_threshold4=0.97;

[behavior_filepaths, temp]=uigetfile('*.mat', 'Chose files to load: ! SAME plane!','MultiSelect','on');
% should find the distance first
[figure_filepaths, fig_temp]=uigetfile('*.fig', 'Chose figure to load:','MultiSelect','on');
f_ROI=openfig(figure_filepaths); 
movegui(f_ROI,'northwest')
for f=1:size(behavior_filepaths,2)
load(behavior_filepaths{f});

leak_id=[];
num_pf_id=[];
pf_count=0;
center_pos=[];
pf_id=[];
pf_simval=[];
  for i=1:size(sig_PFs,2)  
    for j=1:size(sig_PFs,1)
        if ~isempty(sig_PFs{j,i})
             pf_count=pf_count+1;
             pf_id=[pf_id i];
%            sig_PF_M(:,:,pf_count)=sig_PFs{j,i};
             [center_pos(pf_count,1) center_pos(pf_count,2)]=find_center(cell_ROI(:,:,i));  
             
%              figure;
%              imagesc(sig_PFs{j,i}');
%              title(num2str([i j]))
             
        end
    end
  end
  distance= pdist(center_pos);
  dis_M=zeros(size(center_pos,1)-1,size(center_pos,1));
  for i=1:size(center_pos,1)-1
      for j=i+1:size(center_pos,1)
      dis_M(i,j)=distance((i-1)*size(center_pos,1)+j-i*(i+1)/2);
      end
      
      
  end
  
  %% find ROIs that are close to each other
  if ~isempty(dis_M)
   neighbor_id=[]; 
  [neighbor_id(:,1) neighbor_id(:,2)]=find(dis_M<dist_threshold);% & dis_M>0);
  true_neighbor_id=[];
  true_neighbor_id=neighbor_id(find(neighbor_id(:,1)<neighbor_id(:,2)),:);
  

    
  for i=1:size(true_neighbor_id,1)
      
      if i==1
       if pf_id(true_neighbor_id(i,1))~=pf_id(true_neighbor_id(i,2))
      id_pf1=1; id_pf2=1;
     pf1= sig_PFs{id_pf1,pf_id(true_neighbor_id(i,1))}';
     pf2= sig_PFs{id_pf2,pf_id(true_neighbor_id(i,2))}';
     
      else
      id_pf1=1; id_pf2=2;
      pf1= sig_PFs{id_pf1,pf_id(true_neighbor_id(i,1))}';
      pf2= sig_PFs{id_pf2,pf_id(true_neighbor_id(i,2))}';         
          
      end
      
      else
          
      if pf_id(true_neighbor_id(i,1))~=pf_id(true_neighbor_id(i,2)) 
          if pf_id(true_neighbor_id(i-1,1))== pf_id(true_neighbor_id(i,1))&true_neighbor_id(i-1,1)~=true_neighbor_id(i,1)
              if number_of_PFs(true_neighbor_id(i,2))==1
                id_pf1=2; id_pf2=1;
                pf1= sig_PFs{id_pf1,pf_id(true_neighbor_id(i,1))}';
                pf2= sig_PFs{id_pf2,pf_id(true_neighbor_id(i,2))}';
              else
                if pf_id(true_neighbor_id(i,2))==pf_id(true_neighbor_id(i,2)-1)
                id_pf1=2; id_pf2=2;
                pf1= sig_PFs{id_pf1,pf_id(true_neighbor_id(i,1))}';
                pf2= sig_PFs{id_pf2,pf_id(true_neighbor_id(i,2))}'; 
                else 
                id_pf1=2; id_pf2=1;
                pf1= sig_PFs{id_pf1,pf_id(true_neighbor_id(i,1))}';
                pf2= sig_PFs{id_pf2,pf_id(true_neighbor_id(i,2))}';
                end
              end           
          
     
            else
              if pf_id(true_neighbor_id(i,2))== pf_id(true_neighbor_id(i,2)-1)
                id_pf1=1; id_pf2=2;
                pf1= sig_PFs{id_pf1,pf_id(true_neighbor_id(i,1))}';
                pf2= sig_PFs{id_pf2,pf_id(true_neighbor_id(i,2))}';
              else 
                id_pf1=1; id_pf2=1;
                pf1= sig_PFs{id_pf1,pf_id(true_neighbor_id(i,1))}';
                pf2= sig_PFs{id_pf2,pf_id(true_neighbor_id(i,2))}';   
               end
            
          end
        else
        id_pf1=1; id_pf2=2;
        pf1= sig_PFs{id_pf1,pf_id(true_neighbor_id(i,1))}';
        pf2= sig_PFs{id_pf2,pf_id(true_neighbor_id(i,2))}';
          
      end
      end
      
   %ssimval = ssim(pf1,pf2);
   ssimval=find_co(pf_id(true_neighbor_id(i,1)),pf_id(true_neighbor_id(i,2)),id_pf1,id_pf2,sig_PFs);
   pf_simval=[pf_simval ssimval];
   
      id1=true_neighbor_id(i,1);
      id2=true_neighbor_id(i,2);
      pf_id(true_neighbor_id(i,1))
      pf_id(true_neighbor_id(i,2))
      cur_dis=dis_M(id1,id2)
      ssimval
      overlap=0;
%       
%       if cur_dis>dist_step1 & ssimval>sim_threshold1
%           overlap=1;
%           
%       elseif cur_dis>dist_step2 & cur_dis<=dist_step1 & ssimval>sim_threshold2
%           overlap=1;
%           
%       elseif cur_dis>dist_step3 & cur_dis<=dist_step2 & ssimval>sim_threshold3
%           overlap=1; 
%           
%       elseif cur_dis<dist_step3 & ssimval>sim_threshold4
%           overlap=1;
%           
%       end
   if ssimval>sim_threshold1%thres_calculator(base_threshold,dist_threshold,cur_dis)
       overlap=1;
   end
      if overlap==1
          
%           leak_id=[leak_id pf_id(true_neighbor_id(i,2))];
%           num_pf_id=[num_pf_id id_pf2];

         pf1_wn=sig_PFs_with_noise{id_pf1,pf_id(true_neighbor_id(i,1))}';
         pf2_wn=sig_PFs_with_noise{id_pf2,pf_id(true_neighbor_id(i,2))}';
         comp_p=figure;
         subplot(2,2,1); imagesc(pf1);title([num2str(pf_id(true_neighbor_id(i,1))) 'pf=' num2str(id_pf1)]);
         subplot(2,2,2); imagesc(pf1_wn);title([num2str(pf_id(true_neighbor_id(i,1))) 'pf=' num2str(id_pf1)]);colorbar;
         subplot(2,2,3); imagesc(pf2);title([num2str(pf_id(true_neighbor_id(i,2))) 'pf=' num2str(id_pf2)]);    
         subplot(2,2,4); imagesc(pf2_wn);title([num2str(pf_id(true_neighbor_id(i,2))) 'pf=' num2str(id_pf2)]);colorbar;    
          movegui(comp_p,'southwest')
          

          figure;
          imagesc(cell_ROI(:,:,pf_id(true_neighbor_id(i,1)))-cell_ROI(:,:,pf_id(true_neighbor_id(i,2))));
          varified_result = input(['is there leakage between the 2 ROI? 1 for yes 0 for no>>']);
            figure(f_ROI);
            figure(comp_p);
          
          if varified_result==1
          chose_remove = input(['which ROI you wants to remove? 1 or 2>>']);
          if chose_remove~=1 & chose_remove~=2
              chose_remove = input(['You seem putted a wrong input.Which ROI you wants to remove? 1 or 2>>']);
          end
          leak_id=[leak_id pf_id(true_neighbor_id(i,chose_remove))];
          if chose_remove==1
          num_pf_id=[num_pf_id id_pf1];
          else
          num_pf_id=[num_pf_id id_pf2];
          end
%         figure;
%         subplot(2,1,1); imagesc(pf1);title(num2str(pf_id(true_neighbor_id(i,1))));
%         subplot(2,1,2); imagesc(pf2);title(num2str(pf_id(true_neighbor_id(i,2))));    
          varified_result=[]; 
          chose_move=[];
          end
              
      end
  
   
  end  
  for i=1:size(leak_id)
   %% remove leakage inbetween animals
%   if overlap==1   %ssimval>sim_threshold
%       step down the value to pickout possible same pf 

      
%       pf_id(true_neighbor_id(i,1))
%       pf_id(true_neighbor_id(i,2))
%       sim_score=dis_M(true_neighbor_id(i,1),true_neighbor_id(i,2))*ssimval
%       
%       if sim_score<22

    remove_id=leak_id(i);
    remove_pf=num_pf_id(i);
    sig_PFs{remove_pf,remove_id}=[];  
    sig_PFs_with_noise{remove_pf,remove_id}=[]; 
    number_of_PFs(remove_id)=number_of_PFs(pf_id(true_neighbor_id(i,2)))-1; 
    PF_end_bins(remove_pf,remove_id)=nan;
    PF_start_bins(remove_pf,remove_id)=nan;
    PF_PVALS(remove_pf,remove_id)=nan;
    PF_rate(remove_pf,remove_id)=nan;
    PF_width(remove_pf,remove_id)=nan;
    
%     figure;
%     subplot(2,1,1); imagesc(pf1);title(num2str(pf_id(true_neighbor_id(i,1))));
%     subplot(2,1,2); imagesc(pf2);title(num2str(pf_id(true_neighbor_id(i,2))));
    
    
%       end
%   end
      
%   end
    
  end
  end
 save([behavior_filepaths{f}(1:end-4) '_leakageremoved'],'sig_PFs','sig_PFs_with_noise','PF_width','PF_rate','PF_start_bins','PF_PVALS','PF_end_bins','number_of_PFs','mean_trans','label_high','label_low','ybinmax','cell_ROI','leak_id');%,'-append'); 
end

%% usful function
%pdist
