function [x y jump_y jump_x]=linreg_for_bylap_groupedCOM(sig_PFs,pf_id,norm_COMshift,startlap,stepsize,standardlap,stoplap)
y=[];
x=[];
jump_y=[];
jump_x=[];

 for i =1:size(sig_PFs,2)
     for p=1:size(sig_PFs,1)
       cur_sig_PFs=sig_PFs{p,i};
       cur_COMshift=norm_COMshift{p,i};
       use_id=find(cur_COMshift>-10& cur_COMshift<10);  
       %lapnum=size(cur_sig_PFs{1,use_id(1)},2);  
       cur_pf_id=pf_id{p,i};
       for j=use_id
            cur_startlap=startlap{p,i}(j);
            if cur_startlap<20
            lapnum=size(cur_sig_PFs{1,cur_pf_id(j)},2);  
            lapnum=min(lapnum,stoplap+stepsize);
            cur_binmean=cur_sig_PFs{1,cur_pf_id(j)};
            if sum(mean(cur_binmean(:,1:5),2)>0)<51
            cur_COMv=[];
            cur_jump_y=[];
            cur_jump_x=[];
            for n=1:(lapnum-stepsize+1)
                step_binmean=cur_binmean(:,n:n+stepsize-1);
                [meanCOM SP]=meanCOMandSP(step_binmean',1,50,50);
                cur_COMv=[cur_COMv meanCOM];
            if mod(n,stepsize)==1
               cur_jump_y=[cur_jump_y meanCOM];
               cur_jump_x=[cur_jump_x n];
                
                
            end
            end  
            cur_y=cur_COMv-cur_COMv(standardlap);
            cur_x=1:(lapnum-stepsize+1);
            remove_id=isnan(cur_y);
            cur_y(remove_id)=[];                    
            cur_x(remove_id)=[];
            y=[y cur_y];
            x=[x cur_x];
            
            cur_jump_y=cur_jump_y-cur_COMv(standardlap+1);
            jump_remove_id=isnan(cur_jump_y);
            cur_jump_y(jump_remove_id)=[];
            cur_jump_x(jump_remove_id)=[];
            jump_y=[jump_y cur_jump_y];
            jump_x=[jump_x cur_jump_x];
            end
            end    
        end
        
       end
         
     end
    
 end



%  model=fitlm(x,y);
%  
% figure;plot(x,y,'o')
%  range=1:1:max(x);
% [count id]=histc(x,range);
% mean_y=[];
% sem=[];
% for i=range
%     cur_y=y(id==i);
%     mean_y=[mean_y mean(y(id==i))];
%     sem= [sem std(cur_y)/sqrt(length(cur_y))];
% end
% 
% 
% %mdl=fitlm(x,y);
% hold on;
% a=model.('Coefficients').('Estimate')(2);
% plot(x,model.('Coefficients').('Estimate')(2).*x+ model.('Coefficients').('Estimate')(1));
% title(['slope=' num2str(a)]);
% ylim([-50 50]);
% 
% 
% % figure;
% % %plot(range,mean_y,'o')
% % errorbar(range,mean_y,sem,'o')
% % ylim([-10 10])
% end