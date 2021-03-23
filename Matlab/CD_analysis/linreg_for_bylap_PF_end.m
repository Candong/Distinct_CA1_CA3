function [x y x_M y_M]=linreg_for_bylap_PF_end(f_sig_PFs,pf_id,stepsize,standardlap,onsetlap,stoplap,find_start)
% calculate windowed COM of PF acroos all animals 
% x:lap number y: pf_shift
y=[];
x=[];
count=0;
window=6;
threshold=3;
x_M=[];
y_M=[];
stoplap=stoplap+1;
 for i =1:size(f_sig_PFs,2)
     for p=1:size(f_sig_PFs,1)
       cur_sig_PFs=f_sig_PFs{p,i};
       cur_pf_id=pf_id{p,i};
       for j=cur_pf_id'
           binmean=cur_sig_PFs{1,j};
           %cur_startlap=startlap{p,i}(j);
           
           if ~isclipped(binmean)
                if onsetlap~=0
                    cur_startlap=find_delaylap(binmean,1,size(binmean,1),window,threshold,size(binmean,2));
                    binmean=binmean(:,cur_startlap:end);
                end
                %skewness_lap=skewness(binmean);
                lapnum=size(binmean,2);
                end_lap=find_pfend(binmean,stepsize,find_start);
                
                if lapnum-stepsize+1>=standardlap
                if ~isnan(end_lap(standardlap))
                cur_y=end_lap-end_lap(standardlap);
                cur_x=1:lapnum-stepsize+1;
                if length(cur_y)>=stoplap
                    y_M=[y_M; cur_y(1:stoplap)];
                    x_M=[x_M; cur_x(1:stoplap)];
                else
                    new_cur_y=[cur_y ones(1,stoplap-length(cur_y))*NaN];
                    new_cur_x=[cur_x ones(1,stoplap-length(cur_y))*NaN];
                    y_M=[y_M; new_cur_y(1:stoplap)];
                    x_M=[x_M; new_cur_x(1:stoplap)];
                end
                
                
                remove_id=find(isnan(end_lap));
                cur_y(remove_id)=[];                    
                cur_x(remove_id)=[];
                y=[y cur_y];
                x=[x cur_x];
                
                end
                end
                
           else 
               count=count+1;
           end
%        use_id=find(cur_COMshift>-10& cur_COMshift<10);  

        
       end
         
     end
    %disp(count);
 end
