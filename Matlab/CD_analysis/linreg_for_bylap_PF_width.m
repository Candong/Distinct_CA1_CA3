function [x y x_M y_M PF_onset]=linreg_for_bylap_PF_width(f_sig_PFs,pf_id,stepsize,standardlap,onsetlap,stoplap)
% calculate windowed COM of PF acroos all animals 
% x:lap number y: pf_shift
y=[];
x=[];
count=0;
window=6;
threshold=3;
x_M=[];
y_M=[];
PF_onset=[];
 for i =1:size(f_sig_PFs,2)
     for p=1:size(f_sig_PFs,1)
       cur_sig_PFs=f_sig_PFs{p,i};
       cur_pf_id=pf_id{p,i};
       for j=cur_pf_id'
           binmean=cur_sig_PFs{1,j};
           if onsetlap~=0
           cur_startlap=find_delaylap(binmean,1,size(binmean,1),window,threshold,size(binmean,2));
           binmean=binmean(:,cur_startlap:end);
           end
           lapnum=size(binmean,2);
           if ~isclipped(binmean) & (lapnum-stepsize+1)>=  standardlap
                width_lap=find_pfwidth(binmean,stepsize);
                
                cur_startlap=find_delaylap(binmean,1,size(binmean,1),window,threshold,size(binmean,2));

                if mean(width_lap(width_lap~=0))~=0
                %cur_y=width_lap-width_lap(standardlap);
                cur_y=width_lap;%/mean(width_lap(width_lap~=0));
                cur_x=1:(lapnum-stepsize+1);
                PF_onset=[PF_onset cur_startlap];
                matrix_y=cur_y;
                matrix_y(width_lap==0)=NaN;
                if length(cur_y)>=stoplap
%                     matrix_y=cur_y;
%                     matrix_y(width_lap==0)=NaN;
                    y_M=[y_M; matrix_y(1:stoplap)];
                    x_M=[x_M; cur_x(1:stoplap)];
                else
                    new_cur_y=[matrix_y ones(1,stoplap-length(cur_y))*NaN];
                    new_cur_x=[cur_x ones(1,stoplap-length(cur_y))*NaN];
                    y_M=[y_M; new_cur_y(1:stoplap)];
                    x_M=[x_M; new_cur_x(1:stoplap)];
                end
                remove_id=width_lap==0;
                cur_y(remove_id)=[];                    
                cur_x(remove_id)=[];
                y=[y cur_y];
                x=[x cur_x];
                end
                count=count+1;
           else 
               %count=count+1;
           end

        
       end
         
     end
    %disp(count);
 end


