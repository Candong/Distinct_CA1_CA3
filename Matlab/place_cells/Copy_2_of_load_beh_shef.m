function [ session ] = load_beh_shef( varargin )
%LOAD_beh GUI for aligning .tif files to .mat files
%
% PARAMETERS
%   tif_frame_counts ([]) - The number of frames in each .tif to be
%       aligned. If left empty, tif_filepaths must/will be assigned and the
%       frame counts will come from the files.
%   tif_filepaths ([]) - The filepaths to each .tif. If tif_frame_counts
%       also empty, calls GUI for the filepaths.
%   gap_tolerance (2) - The % of the frame period that will count as a gap
%       between continuous aquisitions.
%   outfile ([]) - If not empty, loads the data stored in the file and adds
%       the aligned behavior data to it.
%   channels ({'Y_position','X_position','reward','lick'}) - The channels
%       from the .mat file to be calculated. Do not add 'Image_sync'. These
%       must match the spelling and capitalization from Axon exactly.
%   out_modes ({'mean','mean',{'count','t'},{'count','t'}}) - The methods
%       for aligning each channel in "channels" with the frames.
%           mean - The average of values within the frame, no suffix.
%           count - The count of the steps that fall within the frame. The
%               time is estimated in the center of each step, suffix '_c'.
%           t - The time (in seconds) of the center of each step from the
%               begining of the first aquisition, suffix '_t'.
%
% OUTPUT
%   session - The session object. If passed outfile, then session contains
%       whatever data was in the original file. Otherwise, it contains a
%       field .mat_data, which contains the fields defined by 'channels'
%       and 'out_modes' and the field 't', which contains the mid-frame
%       time (in seconds) for each assigned frame from the beginning of
%       aquisition.
%
% Jason R. Climer, PhD (jason.r.climer@gmail.com) 18 May, 2017

WINDOW = 10;

%% PARSE INPUT
tif_frame_counts = [];
tif_filepaths = [];
beh_filepaths = [];
gap_tolerance = [];
outfile = [];
columns = [];
out_columns = [];
channels = [];
autoApply = [];

%save_name=['chw3']% put the name you want for the tracked data
% Parameters
ip = inputParser;
ip.addParameter('tif_frame_counts',[]);% The counts of frames to be aligned in each video.
ip.addParameter('tif_filepaths',[]);% The path to the .tif files
ip.addParameter('beh_filepaths',[]);% The path to the .mat files
ip.addParameter('gap_tolerance',2)
ip.addParameter('autoApply',false);
ip.addParameter('outfile',[]);
ip.addParameter('channels',{'Y_pos','reward','airpuff','lick'}); % CAN ANNOTATE might want to change the name here
ip.addParameter('out_modes',{'mean','mean',{'count','t'},{'count','t'}});
ip.parse(varargin{:})
for j=fields(ip.Results)'
    eval([j{1} '=ip.Results.' j{1} ';']);
end

for i=1:numel(out_modes)
    if ~iscell(out_modes{i})
        out_modes{i} = {out_modes{i}};
    end
end

% Handle outfile
if ~isempty(outfile)&&exist(outfile,'file')
    load(outfile);
else
    session = struct;
end

channels = [{'Image_sync'},channels];

% Handle beh_filepaths
if isempty(beh_filepaths)% Not given
    [beh_filepaths, temp]=uigetfile('*.mat', 'Chose behavior files to load:','MultiSelect','on');
     if ~iscell(beh_filepaths), beh_filepaths = {beh_filepaths}; end
    beh_filepaths = cellfun(@(x)[temp x],beh_filepaths,'UniformOutput',false);
    beh_filepaths = sort(beh_filepaths);% Ensure filepaths in lexographic order
    clear temp;
end
 if ~iscell(beh_filepaths), beh_filepaths = {beh_filepaths}; end

% Handle tif_frame_counts
if isempty(tif_frame_counts)% Need to calculate from .tifs
    if isempty(tif_filepaths)% Load .tifs if necessary
        [tif_filepaths, temp]=uigetfile('*.tif', 'Chose image files to load:','MultiSelect','on');
        if ~iscell(tif_filepaths), tif_filepaths = {tif_filepaths}; end
        tif_filepaths = cellfun(@(x)[temp x],tif_filepaths,'UniformOutput',false);
    end
    if ~iscell(tif_filepaths), tif_filepaths = {tif_filepaths}; end
    [~,tif_frame_counts] = load_tiffs_fast(tif_filepaths,'end_ind',1);% Count the frames
    clear temp;
end

%% Load .mat files % CAN ANNOTATE probablly could dirrectly start from here % we don't need the x position
% '.val' is because the way we save the data
%[raw_data,si] = cellfun(@(x)behaviorload(x,'channels',channels),beh_filepaths,'UniformOutput',false);
% for i=1:length(beh_filepaths)
% load_path=beh_filepaths{i};
% whole_data=load(load_path);
% rawdata{1}=whole_data.val.G;
% rawdata{2}=whole_data.val.A;
% rawdata{3}=whole_data.val.B;
% rawdata{4}=whole_data.val.C;
% rawdata{5}=whole_data.val.D;
% raw_data{i}=cell2mat(rawdata);
% 
% end
% 
% 
% si={whole_data.val.Tinterval*1000000};%sample interval in usec;

for i=1:length(beh_filepaths)
load_path=beh_filepaths{i};
whole_data=load(load_path);
rawdata{1}=whole_data.G;
rawdata{2}=whole_data.A;
rawdata{3}=whole_data.B;
rawdata{4}=whole_data.C;
rawdata{5}=whole_data.D;
raw_data{i}=cell2mat(rawdata);

end


si={whole_data.Tinterval*1000000};%sample interval in usec;


% keyboard
data = struct;
data = repmat(data,[numel(beh_filepaths) 1]);
for i=1:numel(channels)
    temp = cellfun(@(x)x(:,i),raw_data,'UniformOutput',false);
    [data(:).(channels{i})] = deal(temp{:});
end

%% Initial calculations from raw data

% Add fields to data
fld_placeholder=cell(size(data));
[data(:).sample_period_t] = deal(fld_placeholder{:});% Period between frames (in .mat samples)

[data(:).on_t] = deal(fld_placeholder{:});% Onset time (in .mat samples) for each frame
[data(:).off_t] = deal(fld_placeholder{:});% Offset time (in .mat samples) for each frame
[data(:).mid_t] = deal(fld_placeholder{:});% Frame midpoint (in .mat samples) for each frame

for i=1:numel(data)% For each .mat file
    data(i).on_t = find(diff(data(i).Image_sync>2)==1)+1;% Find the frame onsets
    data(i).off_t = find(diff(data(i).Image_sync>2)==-1)+1;% Find the frame offsets
    
    if data(i).Image_sync(1)>2% Handle first bad frame
        data(i).off_t = data(i).off_t(2:end);
    end
    if data(i).Image_sync(end)>2% Handle last bad frame
        data(i).on_t = data(i).on_t(1:end-1);
    end
    data(i).mid_t = (data(i).on_t+data(i).off_t)/2;% Find the frame midpoints
    
    data(i).sample_period_t = quantile(diff(data(i).mid_t),0.1); % CAN ANNOTATE don't quite understand probabaly the first diff??
end

%%
beh_count = NaN(sum(cellfun(@numel,{data.Image_sync})),1);

% on_t = cat(1,data.on_t);
% off_t = cat(1,data.off_t);
% mid_t = cat(1,data.mid_t);

on_T = NaN(sum(cellfun(@numel,{data.on_t})),1);
off_T = on_T;
mid_T = on_T;

for i=1:numel(data)
    beh_count(sum(cellfun(@numel,{data(1:i-1).Image_sync}))+1:sum(cellfun(@numel,{data(1:i).Image_sync})))=i;
    on_T(sum(cellfun(@numel,{data(1:i-1).on_t}))+1:sum(cellfun(@numel,{data(1:i).on_t})))=...
        data(i).on_t+...
        sum(cellfun(@numel,{data(1:i-1).Image_sync}));
    off_T(sum(cellfun(@numel,{data(1:i-1).on_t}))+1:sum(cellfun(@numel,{data(1:i).on_t})))=...
        data(i).off_t+...
        sum(cellfun(@numel,{data(1:i-1).Image_sync}));
    mid_T(sum(cellfun(@numel,{data(1:i-1).on_t}))+1:sum(cellfun(@numel,{data(1:i).on_t})))=...
        data(i).mid_t+...
        sum(cellfun(@numel,{data(1:i-1).Image_sync}));
end

%%
blocks = [0;find(diff(mid_T)>mean([data.sample_period_t])*(1+gap_tolerance));numel(mid_T)];
blocks = [blocks(1:end-1)+1 blocks(2:end)];
% plot(1:numel(mid_T),mid_T,'b-',blocks(:,1),5e4*ones(size(blocks(:,1))),'gx',blocks(:,2),5e4*ones(size(blocks(:,1))),'rx')
block_duration = blocks(:,2)-blocks(:,1)+1;
%% Now assign the frames
frame_i = NaN(sum(tif_frame_counts),1);% The index (in detected frames from the start) for each image frame

if numel(frame_i)>=numel(mid_T)
    frame_i(:) = 1:numel(frame_i);
else
    for i=1:numel(tif_frame_counts)
        % If this .tif fits perfectly into an EMPTY block that isn't too
        %   late, put it there.
        % If not, put it into the EMPTY block that isn't too late
        %   that matches the best.
        % If there are no EMPTY blocks, put it at the first available frame
        empty_blocks = cumsum(block_duration);
        empty_blocks = [[1;empty_blocks(1:end-1)+1] empty_blocks];
        empty_blocks = empty_blocks(empty_blocks(:,1)>nansum(nanmax(frame_i))&blocks(end)-empty_blocks(:,1)>sum(tif_frame_counts(i+1:end)),:);
        
        if ~isempty(empty_blocks)
            [~,j] = min(abs(diff(empty_blocks,[],2)+1-tif_frame_counts(i)));
            frame_i(sum(tif_frame_counts(1:i-1))+1:sum(tif_frame_counts(1:i)))=empty_blocks(j,1)-1+(1:tif_frame_counts(i));
        else
            frame_i(sum(tif_frame_counts(1:i-1))+1:sum(tif_frame_counts(1:i)))=nanmax(frame_i)+(1:tif_frame_counts(i));
        end
    end
end
%%
Image_sync = cat(1,data.Image_sync)>2;
beh_file_dur = arrayfun(@(x)numel(x.Image_sync),data);

selected_tif = 1;
dropped_frames = zeros(2,numel(tif_frame_counts));

if ~autoApply
%%
% close all;clc;
figh = figure('units','normalized','Position',[0    0.0370    1.0000    0.8917]);
plotax = axes('Position',[0.1 0.3 0.8 0.6]);

uicontrol(figh,'Style','text','FontSize',13 ...
    ,'units','normalized','String','Start'...
    ,'Position',[0.1 0.125 0.05 0.05]);
uicontrol(figh,'Style','text','FontSize',13 ...
    ,'units','normalized','String','End'...
    ,'Position',[0.15 0.125 0.05 0.05]);

uicontrol(figh,'Style', 'pushbutton', 'String', 'Next .tif [W]','FontSize',16 ...
    ,'units','normalized','Position', [0.1 0.225 0.1 0.05],'UserData',0 ...
    ,'Callback', @(source,event)choosetif(1));
uicontrol(figh,'Style', 'pushbutton', 'String', 'Prior .tif [Q]','FontSize',16 ...
    ,'units','normalized','Position', [0.1 0.175 0.1 0.05],'UserData',0 ...
    ,'Callback', @(source,event)choosetif(-1));

begin_drop_frames = uicontrol(figh,'Style','edit','FontSize',16 ...
    ,'units','normalized','String','0'...
    ,'Position',[0.1 0.1 0.05 0.05]...
    ,'Callback', @(source,event)update_dropped_frames);
end_drop_frames = uicontrol(figh,'Style','edit','FontSize',16 ...
    ,'units','normalized','String','0'...
    ,'Position',[0.15 0.1 0.05 0.05]...
    ,'Callback', @(source,event)update_dropped_frames);

uicontrol(figh,'Style', 'pushbutton', 'String', 'Next block [S]','FontSize',16 ...
    ,'units','normalized','Position', [0.205 0.225 0.1 0.05],'UserData',0 ...
    ,'Callback', @(source,event)mv_block_up);
uicontrol(figh,'Style', 'pushbutton', 'String', 'Last block [A]','FontSize',16 ...
    ,'units','normalized','Position', [0.205 0.175 0.1 0.05],'UserData',0 ...
    ,'Callback', @(source,event)mv_block_dn);
uicontrol(figh,'Style', 'pushbutton', 'String', 'Next frame [X]','FontSize',16 ...
    ,'units','normalized','Position', [0.205 0.125 0.1 0.05],'UserData',0 ...
    ,'Callback', @(source,event)mv_frame_up);
uicontrol(figh,'Style', 'pushbutton', 'String', 'Last frame [Z]','FontSize',16 ...
    ,'units','normalized','Position', [0.205 0.075 0.1 0.05],'UserData',0 ...
    ,'Callback', @(source,event)mv_frame_dn);

uicontrol(figh,'Style', 'pushbutton', 'String', 'Zoom extents','FontSize',16 ...
    ,'units','normalized','Position', [0.31 0.225 0.1 0.05],'UserData',0 ...
    ,'Callback', @(source,event)xlim(mean([si{:}])*1e-7*[1 numel(Image_sync)]));
uicontrol(figh,'Style', 'pushbutton', 'String', 'Zoom .tif','FontSize',16 ...
    ,'units','normalized','Position', [0.31 0.175 0.1 0.05],'UserData',0 ...
    ,'Callback', @zm_tif);
uicontrol(figh,'Style', 'pushbutton', 'String', 'Zoom start','FontSize',16 ...
    ,'units','normalized','Position', [0.31 0.125 0.1 0.05],'UserData',0 ...
    ,'Callback', @zm_start);
uicontrol(figh,'Style', 'pushbutton', 'String', 'Zoom end','FontSize',16 ...
    ,'units','normalized','Position', [0.31 0.075 0.1 0.05],'UserData',0 ...
    ,'Callback', @zm_end);

uicontrol(figh,'Style', 'pushbutton', 'String', 'Apply','FontSize',16 ...
    ,'units','normalized','Position', [0.75 0.150 0.1 0.05],'UserData',0 ...
    ,'Callback', @(source,event)apply_fun);

%%
tif_frame_counts=[38995 39000];
drawfig();
xlim(mean([si{:}])*1e-7*[1 numel(Image_sync)])

uiwait(figh);
else
   apply_fun; 
end

    function [] = update_dropped_frames
        if all(isstrprop(get(begin_drop_frames,'String'),'digit'))
            if str2double(get(begin_drop_frames,'String'))>dropped_frames(1,selected_tif)% Dropping more frames from start
                dropped_frames(1,selected_tif) = min(str2double(get(begin_drop_frames,'String')),tif_frame_counts(selected_tif)-dropped_frames(2,selected_tif)-1);
                frame_i(sum(tif_frame_counts(1:selected_tif-1))+(1:dropped_frames(1,selected_tif))) = NaN;
            elseif str2double(get(begin_drop_frames,'String'))<dropped_frames(1,selected_tif)% Adding frames to start
                % Find the earliest available frame before the current
                % beginning of this .tif
                target_start_frame = nanmax([0;frame_i(1:sum(tif_frame_counts(1:selected_tif-1))+1+dropped_frames(1,selected_tif)-1)])+1;
                % Find the frame that would match the requested number of added
                % frames or the earliest available frame, whichever is later
                target_start_frame = max(target_start_frame,...
                    frame_i(sum(tif_frame_counts(1:selected_tif-1))+1+dropped_frames(1,selected_tif))-dropped_frames(1,selected_tif)+str2double(get(begin_drop_frames,'String'))); % Requested frame
                % Translate into new dropped frame
                dropped_frames_temp=dropped_frames(1,selected_tif)-frame_i(sum(tif_frame_counts(1:selected_tif-1))+1+dropped_frames(1,selected_tif))+target_start_frame;
                % Assign new indices to frames
                frame_i(sum(tif_frame_counts(1:selected_tif-1))+1+dropped_frames(1,selected_tif)-(dropped_frames(1,selected_tif)-dropped_frames_temp:-1:1))=...
                    frame_i(sum(tif_frame_counts(1:selected_tif-1))+1+dropped_frames(1,selected_tif))...% current start frame
                    -(dropped_frames(1,selected_tif)-dropped_frames_temp:-1:1);% Added frames
                
                % Assign dropped_frames
                dropped_frames(1,selected_tif) = dropped_frames_temp;
            end
        end
        
        if all(isstrprop(get(end_drop_frames,'String'),'digit'))
            if str2double(get(end_drop_frames,'String'))>dropped_frames(2,selected_tif)% Dropping more frames from end
                dropped_frames(2,selected_tif) = min(str2double(get(end_drop_frames,'String')),tif_frame_counts(selected_tif)-dropped_frames(1,selected_tif)-1);
                frame_i(sum(tif_frame_counts(1:selected_tif))+1-(1:dropped_frames(2,selected_tif))) = NaN;
            elseif str2double(get(end_drop_frames,'String'))<dropped_frames(2,selected_tif)% Adding frames to end
                % Find the latest available frame after the current end of
                % this .tif
                target_frame = nanmin([blocks(end)+1;frame_i(sum(tif_frame_counts(1:selected_tif))-dropped_frames(2,selected_tif)+1:end)-1]);
                % Find the frame that would match the requested number of added
                % frames or the earliest available frame, whichever is
                % earlier
                target_frame = min(target_frame,frame_i(sum(tif_frame_counts(1:selected_tif))-dropped_frames(2,selected_tif))+dropped_frames(2,selected_tif)-str2double(get(end_drop_frames,'String')));
                % Translate into new dropped frame
                dropped_frames_temp = frame_i(sum(tif_frame_counts(1:selected_tif))-dropped_frames(2,selected_tif))+dropped_frames(2,selected_tif)-target_frame;
                % Assign new indicies to frames
                frame_i(sum(tif_frame_counts(1:selected_tif))+(-dropped_frames(2,selected_tif):-dropped_frames_temp)+1)=...
                    frame_i(sum(tif_frame_counts(1:selected_tif))-dropped_frames(2,selected_tif))+(1:dropped_frames(2,selected_tif)-dropped_frames_temp+1);
                % Assign dropped_frames
                dropped_frames(2,selected_tif) = dropped_frames_temp;
            end
            
            
        end
        % Update dropped frames fields
        choosetif(0);
    end

    function drawfig
        axes(plotax)
        cla;
        hold on;
        plot(((1:numel(Image_sync)))*mean([si{:}])*1e-7,Image_sync,'Color',[1 1 1]*0.7);
        xlabel('Aquisition time (sec)');
        ylim([-0.5 1.2]);
        for i=1:numel(data)
            if i<numel(data)
                plot(sum(beh_file_dur(1:i))*[1 1]*mean([si{:}])*1e-7,[0 1],'--','Color',[1 1 1]*0.7);
            end
            text(mean([si{:}])*1e-7*(sum(beh_file_dur(1:i-1))+beh_file_dur(i)/2),-0.05,sprintf('%i/%i',i,numel(data)),'Color',[1 1 1]*0.7,'HorizontalAlignment','center');
        end
        
        colors = lines(numel(tif_frame_counts));
        A=tif_frame_counts;
        for i=1:numel(tif_frame_counts)
            if i==selected_tif
                plot(mean([si{:}])*1e-7*mid_T(frame_i(sum(tif_frame_counts(1:i-1))+(1+dropped_frames(1,i):tif_frame_counts(i)-dropped_frames(2,i)))),0.5*ones(tif_frame_counts(i)-sum(dropped_frames(:,i)),1),'o','MarkerEdgeColor','k','MarkerFaceColor',colors(i,:),'MarkerSize',10);
            else
                plot(mean([si{:}])*1e-7*mid_T(frame_i(sum(tif_frame_counts(1:i-1))+(1+dropped_frames(1,i):tif_frame_counts(i)-dropped_frames(2,i)))),0.5*ones(tif_frame_counts(i)-sum(dropped_frames(:,i)),1),'*','Color',colors(i,:));
            end
            text(mean([si{:}])*1e-7*mean(mid_T(frame_i(sum(tif_frame_counts(1:i-1))+[1+dropped_frames(1,i) tif_frame_counts(i)-dropped_frames(2,i)]))),0.6,sprintf('%i/%i',i,numel(tif_frame_counts)),'Color',colors(i,:),'HorizontalAlignment','center');
            text(mean([si{:}])*1e-7*mean(mid_T(frame_i(sum(tif_frame_counts(1:i-1))+[1+dropped_frames(1,i) tif_frame_counts(i)-dropped_frames(2,i)]))),0.55,...
                sprintf('Frames %i:%i of %i',dropped_frames(1,i)+1,tif_frame_counts(i)-dropped_frames(2,i),tif_frame_counts(i)),'Color',colors(i,:),'HorizontalAlignment','center');
        end
        i = selected_tif;
        scatter(...
            mean([si{:}])*1e-7*mean(mid_T(frame_i(sum(tif_frame_counts(1:i-1))+[1+dropped_frames(1,i) tif_frame_counts(i)-dropped_frames(2,i)]))),0.45,'k^','filled');
        hold off;
        drawnow();
    end

    function [] = zm_tif(varargin)
        xlim(mean([si{:}])*1e-7*(mean(mid_T(frame_i(sum(tif_frame_counts(1:selected_tif-1))+[1+dropped_frames(1,selected_tif) tif_frame_counts(selected_tif)-dropped_frames(2,selected_tif)])))+[-1 1]*0.55*diff(mid_T(frame_i(sum(tif_frame_counts(1:selected_tif-1))+[1+dropped_frames(1,selected_tif) tif_frame_counts(selected_tif)-dropped_frames(2,selected_tif)])))));
    end

    function [] = zm_start(varargin)
        xlim(mean([si{:}])*1e-7*(mid_T(frame_i(sum(tif_frame_counts(1:selected_tif-1))+1+dropped_frames(1,selected_tif)))+[-1 1]*WINDOW*mean([data.sample_period_t])));
    end

    function [] = zm_end(varargin)
        xlim(mean([si{:}])*1e-7*(mid_T(frame_i(sum(tif_frame_counts(1:selected_tif))-dropped_frames(2,selected_tif)))+[-1 1]*WINDOW*mean([data.sample_period_t])));
    end

    function mv_block_up
        % Find the last frame of this session
        last_frame = frame_i(sum(tif_frame_counts(1:selected_tif))-dropped_frames(2,selected_tif));
        % Find the block that we're targeting the last frame into
        if ismember(last_frame,blocks(:,2))
            target_block = min(find(ismember(blocks(:,2),last_frame))+1,size(blocks,1));
        else
            target_block = find(last_frame>=blocks(:,1)&last_frame<=blocks(:,2));
        end
        
        
        % Find the lastest available frame in that block
        new_end_frame = [];
        while isempty(new_end_frame)
            new_end_frame = blocks(target_block,1):blocks(target_block,2);
            new_end_frame = new_end_frame(find(~ismember(new_end_frame,frame_i([1:sum(tif_frame_counts(1:selected_tif-1)) sum(tif_frame_counts(1:selected_tif))+1:end])),1,'last'));
            target_block = target_block - 1;
        end
        if ~isempty(new_end_frame)
            frame_i(sum(tif_frame_counts(1:selected_tif-1))+1+dropped_frames(1,selected_tif):sum(tif_frame_counts(1:selected_tif))-dropped_frames(2,selected_tif)) = ...
                new_end_frame-(tif_frame_counts(selected_tif)-sum(dropped_frames(:,selected_tif)):-1:1)+1;
        end
        drawfig();
    end

    function mv_block_dn
        % Find the first frame of this session
        first_frame = frame_i(sum(tif_frame_counts(1:selected_tif-1))+1+dropped_frames(1,selected_tif));
        % Find the block that we're targeting the first frame into
        if ismember(first_frame,blocks(:,1))
            target_block = max(find(ismember(blocks(:,1),first_frame))-1,1);
        else
            target_block = find(first_frame>=blocks(:,1)&first_frame<=blocks(:,2));
        end
        % Find the earliest available frame in that block
        new_start_frame = [];
        while isempty(new_start_frame)&&target_block<size(blocks,1)
            new_start_frame = blocks(target_block,1):blocks(target_block,2);
            new_start_frame = new_start_frame(find(~ismember(new_start_frame,frame_i([1:sum(tif_frame_counts(1:selected_tif-1)) sum(tif_frame_counts(1:selected_tif))+1:end])),1));
            target_block = target_block+1;
        end
        % Move to new start frame
        if ~isempty(new_start_frame)
            frame_i(sum(tif_frame_counts(1:selected_tif-1))+1+dropped_frames(1,selected_tif):sum(tif_frame_counts(1:selected_tif))-dropped_frames(2,selected_tif)) = ...
                new_start_frame+(1:tif_frame_counts(selected_tif)-sum(dropped_frames(:,selected_tif)))-1;
        end
        drawfig();
    end

    function mv_frame_dn
        % Find if earlier frame is available
        if ~ismember(frame_i(sum(tif_frame_counts(1:selected_tif-1))+1+dropped_frames(1,selected_tif))-1,frame_i)&&frame_i(sum(tif_frame_counts(1:selected_tif-1))+1+dropped_frames(1,selected_tif))-1>0
            frame_i(sum(tif_frame_counts(1:selected_tif-1))+1+dropped_frames(1,selected_tif):sum(tif_frame_counts(1:selected_tif))-dropped_frames(2,selected_tif))=...
                frame_i(sum(tif_frame_counts(1:selected_tif-1))+1+dropped_frames(1,selected_tif):sum(tif_frame_counts(1:selected_tif))-dropped_frames(2,selected_tif))-1;
            drawfig();
        end
    end

    function mv_frame_up
        % Find if earlier frame is available
        if ~ismember(frame_i(sum(tif_frame_counts(1:selected_tif))-dropped_frames(2,selected_tif))+1,frame_i)&&frame_i(sum(tif_frame_counts(1:selected_tif))-dropped_frames(2,selected_tif))+1<blocks(end)
            frame_i(sum(tif_frame_counts(1:selected_tif-1))+1+dropped_frames(1,selected_tif):sum(tif_frame_counts(1:selected_tif))-dropped_frames(2,selected_tif))=...
                frame_i(sum(tif_frame_counts(1:selected_tif-1))+1+dropped_frames(1,selected_tif):sum(tif_frame_counts(1:selected_tif))-dropped_frames(2,selected_tif))+1;
            drawfig();
        end
    end


    function [] = choosetif(in)
        selected_tif = mod(selected_tif+in-1,numel(tif_frame_counts))+1;
        set(begin_drop_frames,'Callback',[]);
        set(end_drop_frames,'Callback',[]);
        set(begin_drop_frames,'String',num2str(dropped_frames(1,selected_tif)));
        set(end_drop_frames,'String',num2str(dropped_frames(2,selected_tif)));
        set(begin_drop_frames,'Callback',@(source,event)update_dropped_frames);
        set(end_drop_frames,'Callback',@(source,event)update_dropped_frames);
        drawfig();
    end

    function apply_fun
        for i=2:numel(channels)
            for j=1:numel(out_modes{i-1})
                switch out_modes{i-1}{j}
                    case 'mean'
                        temp = cat(1,data.(channels{i}));
                        session.beh_data.(channels{i}) = arrayfun(@(a,b)mean(temp(a:b)),on_T(frame_i(~isnan(frame_i))),off_T(frame_i(~isnan(frame_i))));
                    case {'count','t'}
                        temp = cat(1,data.(channels{i}));
                        temp = temp>mean(minmax(temp'));
                        on_times = find(diff(temp)==1)+1;
                        off_times = find(diff(temp)==-1)+1;
                        if temp(1), off_times=off_times(2:end);end
                        if temp(end), on_times=on_times(1:end-1);end
                        t = mean([on_times off_times],2);% WARNING: THIS IS WRONG
                        a=[];b=[];
                        switch out_modes{i-1}{j}
                            case 'count'
                                session.beh_data.([channels{i} '_c'])=arrayfun(@(a,b)sum(t>=a&t<b),on_T(frame_i(~isnan(frame_i))),off_T(frame_i(~isnan(frame_i))));
                            case 't'
                                session.beh_data.([channels{i} '_t'])=t*mean([si{:}])*1e-7;
                        end
                    otherwise
                        warning('out_mode "%s" for channel "%s" not supported. Skipping...',out_modes{i-1}{j},channels{i});
                end
            end
        end
        session.beh_data.fr = 1/(median(diff(mid_T))*mean([si{:}])*1e-7);
        session.beh_data.t = ((1:sum(~isnan(frame_i)))-1)/session.beh_data.fr;
        if ~isempty(outfile)
            save(outfile,'session');
        else
            save([beh_filepaths{1}(1:find(beh_filepaths{1}=='_',1,'last')+3) '_ds.mat'],'session');
        end
        try
        close(figh);
        catch err;
        end
    end
%save([save_name '.mat'],'session')
end

