clear all; %close all;

[filepaths, temp]=sort(uigetfile('*.mat', 'Chose deltaCOM files to load:','MultiSelect','on'));
usecolor1=[0/255.0,113/255.0,188/255.0]; %blue
usecolor2=[1,0,0];%'red'
usecolor3=[102/255.0, 45/255.0, 145/255.0]; %purple
usecolor4=[34/255.0,181/255.0,115/255.0]; % green
usecolor5=[241/255.0,90/255.0,36/255.0]; %orange

% for f=1:size(filepaths,2)
%  if contains(filepaths{f},'_fnovel')
%       if contains(filepaths{f},'CA1_')
%         load(filepaths{f}); 
%         CA1_f_x=f_x;
%         CA1_nday1_x=n_x;
%         CA1_f_y=f_y;
%         CA1_nday1_y=n_y;
%         CA1_f_x_M=f_x_M;
%         CA1_f_y_M=f_y_M;
%         CA1_nday1_x_M=n_x_M;
%         CA1_nday1_y_M=n_y_M;
%         
%       elseif contains(filepaths{f},'CA3_')
%         load(filepaths{f}); 
%         CA3_f_x=f_x;
%         CA3_nday1_x=n_x;
%         CA3_f_y=f_y;
%         CA3_nday1_y=n_y;
%         CA3_f_x_M=f_x_M;
%         CA3_f_y_M=f_y_M;
%         CA3_nday1_x_M=n_x_M;
%         CA3_nday1_y_M=n_y_M;        
%      
%       end
%       
%  elseif contains(filepaths{f},'_novel2day')
%      
%       if contains(filepaths{f},'CA1_')
%         load(filepaths{f}); 
%         CA1_nday2_x=n_x;
%         CA1_nday2_y=n_y;
%         CA1_nday2_x_M=n_x_M;
%         CA1_nday2_y_M=n_y_M;
%       elseif contains(filepaths{f},'CA3_')
%         load(filepaths{f}); 
%         CA3_nday2_x=n_x;
%         CA3_nday2_y=n_y;
%         CA3_nday2_x_M=n_x_M;
%         CA3_nday2_y_M=n_y_M;
%         
%      
%       end
%      
%   
%  end
% end

%%
stoplap=25;
for f=1:size(filepaths,2)
 if contains(filepaths{f},'_fnovel')
      if contains(filepaths{f},'CA1_')
        load(filepaths{f}); 
        CA1_f_x=f_x;
        CA1_nday1_x=n_x;
        CA1_f_x(f_x>stoplap)=[];
        CA1_nday1_x(n_x>stoplap)=[];
        CA1_f_y=f_y;
        CA1_f_y(f_x>stoplap)=[];
        CA1_nday1_y=n_y;
        CA1_nday1_y(n_x>stoplap)=[];
        CA1_f_x_M=f_x_M;
        CA1_f_y_M=f_y_M;
        CA1_nday1_x_M=n_x_M;
        CA1_nday1_y_M=n_y_M;
      elseif contains(filepaths{f},'CA3_')
        load(filepaths{f}); 
        CA3_f_x=f_x;
        CA3_nday1_x=n_x;
        CA3_f_x(f_x>stoplap)=[];
        CA3_nday1_x(n_x>stoplap)=[];
        CA3_f_y=f_y;
        CA3_f_y(f_x>stoplap)=[];
        CA3_nday1_y=n_y;
        CA3_nday1_y(n_x>stoplap)=[];
        CA3_f_x_M=f_x_M;
        CA3_f_y_M=f_y_M;
        CA3_nday1_x_M=n_x_M;
        CA3_nday1_y_M=n_y_M;
       
     
      end
      
 elseif contains(filepaths{f},'_novel2day')
     
      if contains(filepaths{f},'CA1_')
        load(filepaths{f}); 
        CA1_nday2_x=n_x;
        CA1_nday2_y=n_y;
        CA1_nday2_x(n_x>stoplap)=[];
        CA1_nday2_y(n_x>stoplap)=[];
        CA1_nday2_x_M=n_x_M;
        CA1_nday2_y_M=n_y_M;
        
      elseif contains(filepaths{f},'CA3_')
        load(filepaths{f}); 
        CA3_nday2_x=n_x;
        CA3_nday2_y=n_y;
        CA3_nday2_x(n_x>stoplap)=[];
        CA3_nday2_y(n_x>stoplap)=[];
        CA3_nday2_x_M=n_x_M;
        CA3_nday2_y_M=n_y_M;
        
     
      end
     
  
 end
end
%%
sampletimes = 1000;
alpha=0.05;
[CA1_f_slope CA1_f_intcpt CA1_f_up CA1_f_low CA1_f_mean]=resample_CA1(CA1_f_x_M,CA1_f_y_M,CA3_f_x_M,CA3_f_y_M,sampletimes,alpha);
[CA1_nday1_slope CA1_nday1_intcpt CA1_nday1_up CA1_nday1_low CA1_nday1_mean]=resample_CA1(CA1_nday1_x_M,CA1_nday1_y_M,CA3_nday1_x_M,CA3_nday1_y_M,sampletimes,alpha);
[CA1_nday2_slope CA1_nday2_intcpt CA1_nday2_up CA1_nday2_low CA1_nday2_mean ]=resample_CA1(CA1_nday2_x_M,CA1_nday2_y_M,CA3_nday2_x_M,CA3_nday2_y_M,sampletimes,alpha);



%% Calculate different condition CA3 slope and intercept
% CA1_f_y(CA1_f_x==max(CA1_f_x))=[];
% CA1_f_x(CA1_f_x==max(CA1_f_x))=[];
% mode00=fitlm(CA1_f_x, CA1_f_y);
figure; hold on;
plot(CA1_f_up, '--k','HandleVisibility','off');
plot(CA1_f_low, '--k','HandleVisibility','off');
plot(CA1_f_mean);

CA3_f_y(CA3_f_x==max(CA3_f_x))=[];
CA3_f_x(CA3_f_x==max(CA3_f_x))=[];
model0=fitlm(CA3_f_x, CA3_f_y);
plot(CA3_f_x,model0.('Coefficients').('Estimate')(2).*CA3_f_x+ model0.('Coefficients').('Estimate')(1));
title('CA1 CA3 familiar')
legend({'CA1','CA3'})

%%
figure; hold on;
plot(CA1_nday1_up, '--k','HandleVisibility','off');
plot(CA1_nday1_low, '--k','HandleVisibility','off');
plot(CA1_nday1_mean);

% CA3_nday1_y(CA3_nday1_x==max(CA3_nday1_x))=[];
% CA3_nday1_x(CA3_nday1_x==max(CA3_nday1_x))=[];
model1=fitlm(CA3_nday1_x, CA3_nday1_y);
plot(CA3_nday1_x,model1.('Coefficients').('Estimate')(2).*CA3_nday1_x+ model1.('Coefficients').('Estimate')(1));
title('CA1 CA3 novel day1')
legend({'CA1','CA3'})
%%
figure; hold on;
plot(CA1_nday2_up, '--k','HandleVisibility','off');
plot(CA1_nday2_low, '--k','HandleVisibility','off');
plot(CA1_nday2_mean);

% CA3_nday2_y(CA3_nday2_x==max(CA3_nday2_x))=[];
% CA3_nday2_x(CA3_nday2_x==max(CA3_nday2_x))=[];
model2=fitlm(CA3_nday2_x, CA3_nday2_y);
plot(CA3_nday2_x,model2.('Coefficients').('Estimate')(2).*CA3_nday2_x+ model2.('Coefficients').('Estimate')(1));
title('CA1 CA3 novel day2')
legend({'CA1','CA3'})


%% resampling


%[CA1_f_slope CA1_f_intcpt CA1_f_up CA1_low]=resample_CA1(CA1_f_x_M,CA1_f_y_M,CA3_f_x_M,CA3_f_y_M,sampletimes);
figure;
plot(CA1_f_slope,CA1_f_intcpt,'o','color',usecolor4);
hold on;
plot(model0.('Coefficients').('Estimate')(2), model0.('Coefficients').('Estimate')(1),'*','color',usecolor5);
title('CA1 CA3 familiar')
xlabel(['slope'])
ylabel(['intecept'])
legend({'CA1 downsample','CA3'})
box off


figure;
plot(CA1_nday1_slope,CA1_nday1_intcpt,'o','color',usecolor4);
hold on;
plot(model1.('Coefficients').('Estimate')(2), model1.('Coefficients').('Estimate')(1),'*','color',usecolor5);
title('CA1 CA3 novel day1')
xlabel(['slope'])
ylabel(['intecept'])
legend({'CA1 downsample','CA3'})
box off

figure;
plot(CA1_nday2_slope,CA1_nday2_intcpt,'o','color',usecolor4);
hold on;
plot(model2.('Coefficients').('Estimate')(2), model2.('Coefficients').('Estimate')(1),'*','color',usecolor5);
title('CA1 CA3 novel day2')
xlabel(['slope'])
ylabel(['intecept'])
legend({'CA1 downsample','CA3'})
box off


