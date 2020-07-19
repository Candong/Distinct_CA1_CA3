clear all;
close all;

% [behavior_filepaths, temp]=uigetfile('*.mat', 'Chose behavior files to load:','MultiSelect','on');
% load([temp behavior_filepaths]);

Filenames = {};
Filedirs={};
filenum = 1;
filechk = 0;
    while filechk == 0
        [F,D] = uigetfile('*.mat', 'Select Tiff File.  If no more to select, click "Cancel"');
        if isequal(F, 0)
            filechk = 1;
        else
            %CurFileDir = [D F];
            Filedirs{filenum}=D;
            Filenames{filenum} =F;
            filenum = filenum +1;
        end
        
    end
%%     
max_transient=0;
    for i=1:filenum-1
      load([Filedirs{i} Filenames{i}]);
            
      prompt = ['is ' Filenames{i}(1:12) ' a sparse label data? yes=1,no=0 '  ];
      is_sparse= input(prompt);
      if is_sparse==1
          openfig([Filedirs{i} Filenames{i}(1:end-3) 'fig']);
          prompt2 = ['type in cell body id in []'  ];
          cellbody_id= input(prompt2);
      end
      
      transient=data.Fc3~=0;
      if is_sparse==1
          transient_id=transient(:,cellbody_id);
      else
          transient_id=transient;
      end
      transient_shape=bwlabel(reshape(transient_id,[],1));

      
      pop_size=[];
      for j=1:max(transient_shape)
          pop_size(j)=length(find(transient_shape==j));
                    
      end
      if max_transient<max(pop_size)
          max_transient=max(pop_size);
      end
      group_popsize{i}=pop_size;

        
    end
    
    edges=1:20:max_transient;
    for i=1:filenum-1
          [N,edges] = histcounts(group_popsize{i},edges);
          
        figure;
        histogram(group_popsize{i},edges);
        %bar(edges,N,'histc');      
        title(Filenames{i}(1:12));
        
    end
    
% 
% 
% if nargin < 1 && ~exist('Tiffs','var')
%     filechk = 0;
% 
% 
% elseif nargin < 1 && exist('Tiffs','var')
%     Filenames = Tiffs;
% elseif nargin == 1
%     Filenames = varargin{1};
% else
%     error('Error finding tiffs. Please read script instructions for usable options');
% end

% is sp12_Plain_3 a sparse label data? yes=1,no=0 1
% type in cell body id in [][21 63 6 55 54 53 41 44 43 45 49 50 57 33 15]
% is sp12_Plain_4 a sparse label data? yes=1,no=0 1
% type in cell body id in [][1 4 5 6 7 12 13 14 15 26 18 20 17 28 24 27]
