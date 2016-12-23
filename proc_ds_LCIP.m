function avgdata = proc_ds_VETSA(filename,graphics,nobasecal,basestart,baselength);

% USAGE: avgdata = proc_ds_new (filename,graphics,basestart,baselength);
%
% PURPOSE: to clean-up, average and calculate variables of interest
% for raw handheld pupil data compiled in one Excel spreadsheet
%
% Assumptions: spreadsheet must have "ID", "Task (PLR)", and "1st data pt" columns
% Rows that have "Bad data" or "no data" in Task column will be removed from further analysis.
% Any row that has no data should be properly coded in "Task (PLR)" column.
%
% Input:
%   filename = spreadsheet filename (.xls or .xlsx) in current directory
%   graphics = toggle 0/1 if diagnostic graphs are to be displayed
%   nobasecal = toggle 0/1 if ignoring baseline blink-check and calculation
%   basestart = starting point for calculation of baseline, default 1
%   baselength = number of point to calculate baseline, default 5
% Output:
%   avgdata = structure containing averaged waveforms and all stats
%    Excel spreadsheets filename_cleaned.csv, filename_avgnorm.csv,
%    filename_avg.csv and filename_stats.csv
%
% Before running this function add these folders to Matlab path
% addpath(genpath('S:\Data Back-Up\pspan\analysis\matlab'));
%
% Written by Sanja Kovacevic, Ph.D. on 4/21/2014
% 
% Modified on 5/5/2014 by Sanja Kovacevic, Ph.D.
%
% Adapted for new data with 30Hz sampling rate by Jia Guo, Ph.D. on 11/7/2014
%
% Modified by Jeremy Elman on 6/14/2016:
% Adapted for use on Linux by saving out to .csv with cell2csv.m 


% Before running this function add these folders to Matlab path
%addpath(genpath('C:\Users\skovacevic\Desktop\sk\Pupillometry\Pupil_RO1\analysis\matlab'));  % Rm 306
%addpath(genpath('C:\Users\skovacevic\Desktop\sk\Pupillometry\Pupil_HH\proclightreflex')); % Rm 306
%addpath(genpath('S:\Data Back-Up\pspan\analysis\matlab')); % Rm 223/224

% read in excel file with all subjects' data
%cd C:\Users\skovacevic\Desktop\sk\Pupillometry\Pupil_HH  %Rm 306
%filename = 'Light Reflex_raw pupil data CBSTS BL.xls';
%filename = 'Light Reflex_raw pupil data eFROPS BL.xls';
%filename = 'Light Reflex_raw pupil data FROPS BL.xls';

%cd 'S:\Handheld Pupil Data\CBSTS HH Pupil Data\CBSTS BL HH Pupil Data' %Rm 223/224
%filename = 'Light Reflex_raw pupil data CBSTS BL.xls';

%graphics = 0;

if ~exist('nobasecal','var')
    % switch to turn off baseline calculation
    nobasecal = 0; % calculate baseline by default
end

if ~exist('basestart','var')
    % satrting point of data taken as baseline
    basestart = 1; % starting from the first point by default
end

if ~exist('baselength','var')
    % num of timepoints taken as baseline
    baselength = 5; % 5 samples are taken as baseline by default
end

% num of std dev for determining suspect trials
stdcutoff = 4;

%averaged waveform filename
s = strfind(filename,'.csv');
if nobasecal % no baseline calculation
    pref = 'NoBASE_';
else
    pref = [];
end
filenameavg = sprintf('%s%s_avg.csv',pref,filename(1:s-1));
filenameavgnorm = sprintf('%s%s_avgnorm.csv',pref,filename(1:s-1));
filenameclean = sprintf('%s%s_cleaned.csv',pref,filename(1:s-1));
filenamest = sprintf('%s%s_stats.csv',pref,filename(1:s-1));

alldata = dataset('File', filename, 'Delimiter', ',');

Nrecordings = size(alldata,1); % number of valid recordings

header = alldata.Properties.VarNames(1:8);

IDcol = strmatch('ID',header);
ID = alldata.(IDcol); 
if isa(ID{1},'numeric') % in case ID is number, convert them to strings
    for i=1:Nrecordings
        ID{i} = num2str(ID{i});
    end
end
TRcol = strmatch('Task',header);
TR = alldata.(TRcol);

PupTTcol = find(strcmpi('Task_DS_',header)); %contains tags for data removal
Data1col = find(strcmpi('x1stDataPt', header)); %column with first data point

% Figure out which data this is
if strfind(header{TRcol},'PLR')
    fprintf('Processing Light Reflex data...\n');
    sf = 30; % sampling frequency (Hz)
    ns = 5;
    conds = 1;
    condsamp = {[1:150]};
elseif strfind(header{TRcol},'DS')
    sf = 30;  % sampling frequency (Hz) for new data
    ns = 12;
    DScol = strmatch('BxTrial',header);
    DS = alldata.(DScol);
    conds = [3 6 9];
%     condsamp = {[30:39]; [59:68]; [88:97]}; % relevant samples for each DS
    condsamp = {[sf*3+1:sf*4]; [sf*6+1:sf*7]; [sf*9+1:sf*10]}; % relevant samples for each DS
    fprintf('Processing Digit Span data...\n');
end



% Trials tagged with no data or bad data will be removed
Bad = find( strcmpi('no data',alldata.(PupTTcol)) | strcmpi('bad data', alldata.(PupTTcol)));
Good = setdiff([1:length(ID)],Bad');

% Keeping only good trials
ID = ID(Good);
TR = TR(Good);
if strfind(header{TRcol},'DS')
    DS = DS(Good);
end;

if strfind(header{TRcol},'PLR')
    data.RawPupilTrials = double(alldata(Good,Data1col:Data1col+floor(5*sf)-1)); % only process data recorded in 0-5s span
elseif strfind(header{TRcol},'DS')
    data.RawPupilTrials = double(alldata(Good,Data1col:Data1col+floor(10.2*sf)-1)); % only process data recorded in 0-10.2s span
end

% if size(data.RawPupilTrials,2)~=750
%     warning('Different number of data points found than expected. Check excel file!')
% end
data.PupilTrials = zeros(size(data.RawPupilTrials));
data.BlinkTrials = zeros(size(data.RawPupilTrials));


% Organize data for processing
numpts = size(data.RawPupilTrials,2);
data.RescaleFactor=1/sf;  %sampling frequency
data.TrialSeconds=linspace(0,1/sf*numpts,numpts)';

% Processing each trial separately
for i = 1:length(ID)
    tempdata.RescaleData = data.RawPupilTrials(i,:);
    %data=puppreprocess(data,.3,.2,1.0,10,7); %PSPAN: passes 7pt weight, 5 is default
    %SK 4/28/2014: will have to change these parameters as they were meant to be used for 60 Hz data
    %regThreshold=0.3,slopeThreshold=0.2,extremeThreshold=1,numSamples=10,weightingWidth=7
    %using 38.6/60*10 = 6 samples instead of 10 to capture the same amount
    %of time
    %using 38.6/60*7 = 5 pt weight for smoothing
    if any(isnan(tempdata.RescaleData(1:50)))
        error(sprintf('Check xls file! %s data should be marked properly!',ID{i}));
    end
%     tempdata=puppreprocess(tempdata,.3,.2,1.0,min(floor(sf/60*10),2),min(floor(sf/60*7),2)); %adjusted for slower sampling rate
    tempdata = stublinks_new(tempdata,0,0,[],0.1,sf);
    
    data.PupilTrials(i,:) = tempdata.NoBlinks;
    data.BlinkTrials(i,:) = tempdata.BlinkTimes;
end

if nobasecal==1 % ignore baseline calculation
    data.baseline = zeros(length(ID),numpts);
else
    %Get baseline for each trial
    data.baseline = repmat(mean(data.PupilTrials(:,basestart:basestart+baselength-1),2),1,numpts);
end
% creates matrix of baselines as first baselength samples of each trial
% Get normed data (leave zeros and 10 as is)
data.NormedPupTrials=data.PupilTrials-data.baseline.*(1-(data.PupilTrials==0 | data.PupilTrials==10));


% Get averaged data for each subject
% subjects = unique(ID,'stable');
subjects = unique(ID);

% Will store stats in avgdata structure
avgdata.id = {};
avgdata.ntr = zeros(length(subjects),length(conds));
avgdata.NormedPupData = zeros(length(subjects),numpts,length(conds));
avgdata.CleanPupData = zeros(length(subjects),numpts,length(conds));
avgdata.Time = [0:1:numpts-1]./sf;

for c = 1:length(conds)
    ntr = sprintf('ntr%i',conds(c));
    avgdata.stats.(ntr) = zeros(length(subjects),1);
    AvgAmpBaseline = sprintf('AvgAmpBaseline%i',conds(c));
    avgdata.stats.(AvgAmpBaseline) = zeros(length(subjects),1);
    AvgAmpDS = sprintf('AvgAmpDS%i',conds(c));
    avgdata.stats.(AvgAmpDS) = zeros(length(subjects),1);
    for t = 1:10
        AvgAmpTW = sprintf('AvgAmpTW%i_%i',t,conds(c));
        avgdata.stats.(AvgAmpTW) = zeros(length(subjects),1);
    end;
    for t = 1:10
        AvgNormAmpTW = sprintf('AvgNormAmpTW%i_%i',t,conds(c));
        avgdata.stats.(AvgNormAmpTW) = zeros(length(subjects),1);
    end
end



% detect bad trials, average across good trials, and get data organized for stats for each subject
for s = 1:length(subjects)
    tempdata = [];
    ind = strmatch(subjects(s),ID);
    
    fprintf('Processing %s... \n',subjects{s});
    fprintf('Found %i trials:', length(ind));
    fprintf(' %i',ind); fprintf('\n');
    
    if length(ind) > ns
        warning(sprintf('Subject %s has more than %i trials. Check the xls file!',ID{ind(1)},ns))
    elseif length(ind) < 2
        warning(sprintf('Subject %s has less than 2 trials. This subject will be removed!', ID{ind(1)}))
    end
    
    tempdata.RawPupilTrials = data.RawPupilTrials(ind,:);
    tempdata.PupilTrials = data.PupilTrials(ind,:);
    tempdata.NormedPupTrials = data.NormedPupTrials(ind,:);
    tempdata.BlinkTrials = data.BlinkTrials(ind,:);
    tempdata.baseline = data.baseline(ind,:);
    
    % Detect bad trials
    if length(ind)>1
        % trials that are more than 4 std above the average at any time point are marked suspect
        threshold = repmat((stdcutoff*std(tempdata.NormedPupTrials)+mean(abs(tempdata.NormedPupTrials))),length(ind),1);
        tempdata.Suspect= any(abs(tempdata.NormedPupTrials)>threshold,2);
        
        drops=tempdata.Suspect;
        if graphics
            if ~exist('h1','var')
                h1=figure('Name','Diagnostic plot');
            else figure(h1); clf;
            end;
            subplot(4,1,1); bar(tempdata.Suspect); ylabel('suspect');
            title(subjects(s));
        end
        
        % drop short trials
        if strfind(header{TRcol},'PLR')
            drops=drops+(sum(tempdata.RawPupilTrials(:,1:124)==0,2)>1);
            tempdata.cond = ones(1,length(drops));
            if graphics subplot(4,1,3); bar(sum(tempdata.RawPupilTrials(:,1:124)==0,2)>1); ylabel('short trial'); end
        elseif strfind(header{TRcol},'DS')
            dropsDS = zeros(length(drops),1);
            tempdata.cond = zeros(1,length(drops));
            
            indDS = strmatch('S3',DS(ind));
            tempdata.cond(indDS) = 3;
            % drop trials with >50% blinks in the 0-4s (4*sf samples) range
            % for condition 3
            dropsDS(indDS) = sum(tempdata.BlinkTrials(indDS,1:4*sf),2)> floor(4*sf/2);
%             dropsDS(indDS) = (sum(tempdata.RawPupilTrials(indDS,condsamp{find(conds==3)})==0,2)> floor(sf/2));
            
            indDS = strmatch('S6',DS(ind));
            tempdata.cond(indDS) = 6;
            % drop trials with >50% blinks in the 0-7s (7*sf samples) range
            % for condition 6
            dropsDS(indDS) = sum(tempdata.BlinkTrials(indDS,1:7*sf),2)> floor(7*sf/2);
%             dropsDS(indDS) = (sum(tempdata.RawPupilTrials(indDS,condsamp{find(conds==6)})==0,2)> floor(sf/2));
            
            indDS = strmatch('S9',DS(ind));
            tempdata.cond(indDS) = 9;
            % drop trials with >50% blinks in the 0-10s (10*sf samples) range
            % for condition 9
            dropsDS(indDS) = sum(tempdata.BlinkTrials(indDS,1:10*sf),2)> floor(10*sf/2);
%             dropsDS(indDS) = (sum(tempdata.RawPupilTrials(indDS,condsamp{find(conds==9)})==0,2)> floor(sf/2));
            
            if graphics subplot(4,1,2); bar(dropsDS); ylabel('> 50% blinks in 0-X critical period'); end
            drops = drops+dropsDS;
        end
        
        if nobasecal==0 % baseline calculation required
            % drop blinks around the onset of the stimulus
            drops=drops+(sum(tempdata.BlinkTrials(:,basestart:basestart+baselength-1),2)> floor(baselength/2));
            if graphics subplot(4,1,4); bar(sum(tempdata.BlinkTrials(:,basestart:basestart+baselength-1),2)> floor(baselength/2)); ylabel('blink at onset'); end
        else
            if graphics subplot(4,1,4); title('baseline calculation ignored'); end
        end
        
        tempdata.drops=drops>0;
    else
        tempdata.drops = ones(size(ind)); %if only one trial, it will be dropped
        tempdata.cond = str2num(DS{ind}(2));
        
        if graphics
            figure(h1); clf; subplot(4,1,1); bar(tempdata.drops); ylabel('dropped');
            title(subjects(s));
        end;
    end;
    
    %loop through conditions
    for c = unique(tempdata.cond)
        indDS = find(tempdata.cond==c);
        good = ind(indDS(find(tempdata.drops(indDS)<1)));  %index into all data
        bad = find(tempdata.drops(indDS)); %index into one subject data
        
        
        avgdata.id(s) = ID(ind(1));
        ntr = sprintf('ntr%i',c);
        avgdata.stats.(ntr)(s) = length(good);
        
        if ~isempty(good) && length(good)>= 1 %accepting average for at least 1 trial for DS
            % get averaged waveform for remaining good trials
            avgdata.ntr(s,find(conds==c)) = length(good);
            avgdata.NormedPupData(s,:,find(conds==c)) = mean(data.NormedPupTrials(good,:),1);
            avgdata.CleanPupData(s,:,find(conds==c)) = mean(data.PupilTrials(good,:),1);
            
            % baseline calculation
            AvgAmpBaseline = sprintf('AvgAmpBaseline%i',c);
            if nobasecal==1 % if baseline calculation is ignored
                avgdata.stats.(AvgAmpBaseline)(s) = 0;
            else % if baseline calculation is required
                basedata = data.PupilTrials(good,basestart:basestart+baselength-1);
                v = nonzeros(basedata);
                if ~isempty(v)
                    avgdata.stats.(AvgAmpBaseline)(s) = mean(v);
                end
            end
            
            condata = data.PupilTrials(good,condsamp{find(conds==c)});
            v = nonzeros(condata);
            AvgAmpDS = sprintf('AvgAmpDS%i',c);
            if ~isempty(v)
                avgdata.stats.(AvgAmpDS)(s) = mean(v);
            end
            
            for t = 1:10
                AvgAmpTW = sprintf('AvgAmpTW%i_%i',t,c);
                AvgNormAmpTW = sprintf('AvgNormAmpTW%i_%i',t,c); %Amplitude normalized to baseline
                condata = data.PupilTrials(good,find(avgdata.Time>=t-1 & avgdata.Time<t));
                v = nonzeros(condata);
                if ~isempty(v)
                    avgdata.stats.(AvgAmpTW)(s) = mean(v);
                    avgdata.stats.(AvgNormAmpTW)(s) = avgdata.stats.(AvgAmpTW)(s)-avgdata.stats.(AvgAmpBaseline)(s);
                end;
            end;

            
        end;
    end;
    
    % if graphics is requested plot all trials for this subject and mark
    % trials that will be dropped
    
    if graphics
        
        for icond = 1:3 % 3 conditions ('3/6/9') in DS experiment
            if ~exist('h2','var')
                h2=figure('Name','Trials','Color','w','Position',[1101 60 560 891]);
            else
                figure(h2)
                clf;
            end;
            
            for t = icond*4-3:icond*4 %size(tempdata.NormedPupTrials,1)
                isubplot = mod(t,4);
                if isubplot ==0
                    isubplot = 4;
                end
                subplot(5,1,isubplot) % each condition has 4 trials
                plot(data.TrialSeconds, tempdata.BlinkTrials(t,:),'Color',[255 153 153]/255); hold on;
                plot(data.TrialSeconds, tempdata.RawPupilTrials(t,:), 'g')
                if strfind(header{TRcol},'PLR')
                    titletxt = sprintf('%s Tr%i',ID{ind(t)},t);
                elseif strfind(header{TRcol},'DS')
                    titletxt = sprintf('%s %s',ID{ind(t)},DS{ind(t)});
                end;
                
                if ismember(t,find(tempdata.drops<1))
                    plot(data.TrialSeconds, tempdata.PupilTrials(t,:),'k')
                    title(sprintf('%s Ok',titletxt),'Color','k')
                else
                    plot(data.TrialSeconds, tempdata.PupilTrials(t,:),'r')
                    title(sprintf('%s BAD',titletxt),'Color','r')
                end
                
                if ~isempty(find(tempdata.drops<1))
                    subplot(5,1,5) % show mean for each condition
                    plot(data.TrialSeconds,squeeze(avgdata.CleanPupData(s,:,:)))
                    titletxt = sprintf('%s Average N=',avgdata.id{s});
                    titletxt1 = sprintf('%i ',avgdata.ntr(s,:));
                    title([titletxt titletxt1]);
                end
                
            end
            % pause to allow time to check the figure for each subject,
            % each condition
            if graphics
                pause
            end
        end

    end;
    
    %report dropped trials for this subject
    fprintf('Dropped trials for %s: ',avgdata.id{s});
    if length(find(tempdata.drops<1))==length(ind)
        fprintf ('none \n');
    elseif length(find(tempdata.drops==1))>0
        fprintf ('%i ',find(tempdata.drops==1));
        fprintf('\n');
    end
    
    
end;

% print averaged waveforms to xls file
xlsdata(1,1:3) = {'ID' 'DS' 'NTr'};
xlsdata(1,4:numpts+3)=num2cell(avgdata.Time);
for c = 1: length(conds)
    xlsdata(1+(c-1)*(length(avgdata.id))+1:1+c*(length(avgdata.id)),1) = avgdata.id;
    xlsdata(1+(c-1)*(length(avgdata.id))+1:1+c*(length(avgdata.id)),2) = num2cell(conds(c));
    xlsdata(1+(c-1)*(length(avgdata.id))+1:1+c*(length(avgdata.id)),3) = num2cell(avgdata.ntr(:,c));
    xlsdata(1+(c-1)*(length(avgdata.id))+1:1+c*(length(avgdata.id)),4:numpts+3) = num2cell(avgdata.CleanPupData(:,:,c));
end;
cell2csv(filenameavg,xlsdata);

% print averaged normalized waveforms to xls file
xlsdata(1,1:3) = {'ID' 'DS' 'NTr'};
xlsdata(1,4:numpts+3)=num2cell(avgdata.Time);
for c = 1: length(conds)
xlsdata(1+(c-1)*(length(avgdata.id))+1:1+c*(length(avgdata.id)),1) = avgdata.id;
xlsdata(1+(c-1)*(length(avgdata.id))+1:1+c*(length(avgdata.id)),2) = num2cell(conds(c));
xlsdata(1+(c-1)*(length(avgdata.id))+1:1+c*(length(avgdata.id)),3) = num2cell(avgdata.ntr(:,c));
xlsdata(1+(c-1)*(length(avgdata.id))+1:1+c*(length(avgdata.id)),4:numpts+3) = num2cell(avgdata.NormedPupData(:,:,c));
end;
cell2csv(filenameavgnorm,xlsdata);

% print out cleaned data to xls file
clear xlsdata;
xlsdata(1,1:length(header)) = header;
xlsdata(2:length(ID)+1,1) = ID;
xlsdata(2:length(ID)+1,TRcol) = TR;
if strfind(header{TRcol},'DS')
    xlsdata(2:length(ID)+1,DScol) = DS;
end
xlsdata(2:length(ID)+1,length(header):numpts+length(header)-1) = num2cell(data.PupilTrials);
cell2csv(filenameclean,xlsdata);

% print stats to xls file
clear xlsdata
stats = avgdata.stats;
xlsdata(1,1) = {'ID'};
varnames = fieldnames(stats);
xlsdata(1,2:length(varnames)+1)=varnames;
xlsdata(2:length(avgdata.id)+1,1) = avgdata.id;
for s = 1: length(varnames)
    xlsdata(2:length(avgdata.id)+1,s+1) = num2cell(getfield(stats,varnames{s})');
end
cell2csv(filenamest,xlsdata);


