%2dsegalazyer
clearvars;
close all;
file2.slashtype = '/';
codedir = cd;
cleandirname = strfind(codedir, file2.slashtype);
file2.codeparent = codedir(1:cleandirname(end))

addpath([file2.codeparent 'ExchangeTools' file2.slashtype 'plot3k']);
addpath([file2.codeparent 'ExchangeTools' file2.slashtype 'DrosteEffect-BrewerMap-a77e675'])
addpath([file2.codeparent 'ExchangeTools' file2.slashtype 'raacampbell-notBoxPlot-2fbf98c' file2.slashtype 'code']);
slashtype = '/';


%% User Inputs:
pathlist = {
    
%  '/Users/draina/Desktop/2018_E1_3_otx2Fixed/Nuclear Intensities/FixLive/KO_N2C_000'
'/Users/draina/Desktop/2018_E1_3_otx2Fixed/Nuclear Intensities/FixLive/KO_N2C_005'
'/Users/draina/Desktop/2018_E1_3_otx2Fixed/Nuclear Intensities/FixLive/KO_N2C_100'
'/Users/draina/Desktop/2018_E1_3_otx2Fixed/Nuclear Intensities/FixLive/WT_ESL'


};

pathlist_labels = {
    
%'KO 0ng'
'KO 5ng'
'KO 100ng'
'WT ESL'
};


%Order Channel Labels according to image
ChannelLabel = {
    'Nanog'
    'Crap'
    'Crap'
    'OTX2'
    
    };

ChannelCalcs = {                        %1 - Nuc, 2-Cyt, 3- Whole Cell, 4 - Nuc/Cyt ratio
    [1,0,0,0]
    [0,0,0,0]
    [0,0,0,0]
    [1,0,0,0]
    };
ff= cell2mat(ChannelCalcs);

calclbl = {
    'Nuclear'
    'Cytoplasm'
    'Whole Cell'
    'Nuc-to-Cyt Ratio'
    
    };

boxplot = 1
conscat = 1
chanscat = 1


% Reorder Channels according to how you want to plot them (expects 4
% channels, so just duplicate or set one to garbage channel). Also expects
% Ch1 to be DAPI - doesn't quantify this channel!
ReorderChan = [ 3; 1; 4; 2];

plotxlabel = [ChannelLabel(ReorderChan)];

for qq = 1:length(ReorderChan)
    ChannelCalcsR{qq,1} = ChannelCalcs{ReorderChan(qq)};
end

%% Main


%Reassign ChannelVec to easy to read vars
chCalc = cell2mat(ChannelCalcs(ReorderChan));

ch1 = ReorderChan(1);
ch2 = ReorderChan(2);
ch3 = ReorderChan(3);
ch4 = ReorderChan(4);
restable = table;

for ctr2 = 1:length(pathlist)
    file(ctr2).path    = pathlist{ctr2};
    file(ctr2).name    = file(ctr2).path(find(file(ctr2).path==slashtype, 1, 'last')+1:end);
    file(ctr2).maindir = file(ctr2).path(1:(find(file(ctr2).path==slashtype, 1, 'last')-1));
    
    main_dir = dir(fullfile(file(ctr2).path));
    for mm = 4:length(main_dir)
        if main_dir(mm).isdir
            
            
            
            %Look for *.txt and .csv files
            file2.cpath = [file(ctr2).path slashtype main_dir(mm).name];
            dir_cell_fix = [dir(fullfile(file2.cpath, '*_fix.txt')); dir(fullfile(file2.cpath, '*_fix.csv'))];
            dir_cell_live= dir(fullfile(file2.cpath, '*_live.txt'));
            
            if ~isempty(dir_cell_fix)&&~isempty(dir_cell_live)
                
                %Find Fixed tag
                switch dir_cell_fix.name(end-2:end)
                    case('txt')
                        delimiter = '\t';
                    case('csv')
                        delimiter = 'comma';
                end
                temptablefix      = readtable([file2.cpath slashtype dir_cell_fix.name], 'delimiter', delimiter);
                
                %Find Live tag
                switch dir_cell_live.name(end-2:end)
                    case('txt')
                        delimiter = '\t';
                end
                temptablelive     = readtable([file2.cpath slashtype dir_cell_live.name], 'delimiter', delimiter);
                
                %Delete unwanted columns + rename in Fixed table
                tablefix2 = [temptablefix.Mean(temptablefix.Ch==ch2), temptablefix.Mean(temptablefix.Ch==ch3),temptablefix.TrackID(temptablefix.Ch==ch2)];
                tablefix3 = table(tablefix2(:,1), tablefix2(:,2), tablefix2(:,3), 'VariableNames', {ChannelLabel{ch2} ChannelLabel{ch3} 'fixedTrackID'});
                
                %Remove fixed trackIDs without corresponding live trackIDs
                fixfinder = ismember(temptablefix.TrackID(temptablefix.Ch==1),temptablelive.TrackID);  
                del2 = find(fixfinder<1);
                tablefix3(del2,:) = [];
                tablefix = sortrows(tablefix3, 'fixedTrackID', 'ascend');
                
                
                %Remove live TrackIDs without corresponding fixed TrackIDs
                livefinder = ismember(temptablelive.TrackID,temptablefix.TrackID(temptablefix.Ch==1));
                del = find(livefinder<1);
                tablelive1 = temptablelive;
                tablelive1(del,:) = [];
                tablelive = sortrows(tablelive1, 'TrackID', 'ascend');
                
                %Error Checking
                if height(tablelive)~=height(tablefix)
                    msgbox('Danger Will Robinson!')
                    keyboard
                end
                
                %Concatenate tables
                tableconcat = [tablefix tablelive];
                
                %Add additional calculations here
                addcol1tit = 'NanogbyOtx2';
                addcol1 = table(tableconcat.Nanog./tableconcat.OTX2, 'VariableNames', {addcol1tit});
                tableconcat = [tableconcat addcol1];
                
                %Accumulate results from the loop
                restable = [restable; tableconcat];
            end
        end
    end
    
    %Accumulate results across treatments
    restotal{ctr2} = restable;
    clear restable
    restable = table;  %reinitialize for concat
    
end


if chanscat ==1
    
    %% Per Channel Scatter
    clear chlabel
    plotflag.type = 'SingleScatter';
    plotflag.corrprint =0;
    file2.savename = 'otxpeak';

    xchan = 2;
    ychan = 12;
    lims.x = [0 70];
    lims.y = [0 0.15];
    
    tablab = restotal{1}.Properties.VariableNames;
    calclbl = '';
    xcalcType = '_';
    ycalcType = '_';
    chlabel{1,1} = char(tablab{xchan});
    chlabel{1,2} = char(tablab{ychan});
    chlabel{2,1} = char(xcalcType);
    chlabel{2,2} = char(ycalcType);
    
    resvec = 0;
    
    for nn = 1:length(restotal)
        
        scatx = table2array(restotal{nn}(:,xchan));
        scaty = table2array(restotal{nn}(:,ychan));
        plotflag.singTreat = nn;
        
        IF_ncplot(plotflag, resvec,scatx, scaty, char(pathlist_labels{nn}), chlabel, calclbl, file, file2, slashtype, lims)
    end
    
    
end


if conscat ==1
    clear chlabel scatx scaty
    plotflag.type = 'ConScatter';
    
    xchan = 1;
    ychan = 2;
    file2.savename = 'nanogotx';
    lims.x = [0 70];
    lims.y = [0 70];
    
    TreatmentList = [1 2 3];
    tablab = restotal{1}.Properties.VariableNames;
    xcalcType = '_';
    ycalcType = '_';
    
    chlabel{1,1} = char(tablab{xchan});
    chlabel{1,2} = char(tablab{ychan});
    chlabel{2,1} = char(xcalcType);
    chlabel{2,2} = char(ycalcType);
    
    resvec = 0;
    
    
    for qq = 1:size(TreatmentList,2)
        scatx{qq} = table2array(restotal{qq}(:,xchan));
        scaty{qq} = table2array(restotal{qq}(:,ychan));
    end
    
    plotflag.conTreat = TreatmentList;
    plotflag.conColour = length(file); %For colours to stay constant
    IF_ncplot(plotflag, resvec,scatx, scaty, pathlist_labels, chlabel, calclbl, file, file2, slashtype, lims)
    
    
end



%
% if boxplot==1
%
%     %% Per-Channel Box Plots (i.e. across all treatment groups for one channel)
%     activeChans = cellfun(@(x) sum(x),ChannelCalcsR)>0; %Reads the ChannelCalcs variable to decide which channels need calcs
%     plotflag = [1 0 0];
%     scatx = 0;
%     scaty = 0;
%     for rr = 1:length(activeChans) %Channels
%         if activeChans(rr)>0
%             for cc = 1:length(ChannelCalcsR{rr}) %Calculations
%                 calcType = ChannelCalcsR{rr}(cc);
%                 if calcType>0
%
%                     switch rr
%                         case 2
%                             resvec = resvec_calc2(calcType,:);
%                         case 3
%                             resvec = resvec_calc3(calcType,:);
%                         case 4
%                             resvec = resvec_calc4(calcType,:);
%                     end
%
%                     chlabel = ChannelLabel(ReorderChan(rr));
%                     calclbl2 = char(calclbl{calcType});
%                     IF_ncplot(plotflag, resvec,scatx, scaty, pathlist_labels, chlabel, calclbl2, file, slashtype)
%                 end
%             end
%         end
%     end
% end
%
%
%   plotflag.type = 'boxplot';
%     scatx = 0;
%     scaty = 0;
%     calcType = [1,2,3,4,5];
%             for cc = calcType %Calculations
%                 if cc>0
%                     chlabel = 'null';
%                     resvec = resvec_calc(cc,:);
%                     calclbl2 = char(calcLabels(cc));
%                     IF_ncplot(plotflag, resvec,scatx, scaty, pathlist_labels, chlabel, calclbl2, file, slashtype)
%                 end
%             end
%





