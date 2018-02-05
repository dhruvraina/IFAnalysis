%2dsegalazyer
clearvars;
close all;

addpath('/Users/draina/Documents/MATLAB/dr_progs/MPI_Progs/ExchangeTools/raacampbell-notBoxPlot-2fbf98c/code');
addpath('/Users/draina/Documents/MATLAB/dr_progs/MPI_Progs/ExchangeTools/DrosteEffect-BrewerMap-a77e675')
addpath('/Users/draina/Documents/MATLAB/dr_progs/MPI_Progs/IFAnalysis');

slashtype = '/';


%% User Inputs:
pathlist = {
    % '/Users/draina/Desktop/2017_eKTR_ERK_Analysis/PDO3_Release'
    % '/Users/draina/Desktop/2017_eKTR_ERK_Analysis/PDO3'
    % '/Users/draina/Desktop/2017_eKTR_ERK_Analysis/ESL'
    
%     '/Users/draina/Desktop/2017_E2_2_June_ReleaseKinetics/1minRelease'
%     '/Users/draina/Desktop/2017_E2_2_June_ReleaseKinetics/2minRelease'
%     '/Users/draina/Desktop/2017_E2_2_June_ReleaseKinetics/3_5minRelease'
%     '/Users/draina/Desktop/2017_E2_2_June_ReleaseKinetics/7_5minRelease'
%     '/Users/draina/Desktop/2017_E2_2_June_ReleaseKinetics/15minRelease'
%     '/Users/draina/Desktop/2017_E2_2_June_ReleaseKinetics/30minRelease'
%     '/Users/draina/Desktop/2017_E2_2_June_ReleaseKinetics/ESL'
%     '/Users/draina/Desktop/2017_E2_2_June_ReleaseKinetics/PDO3'
    
    '/Users/draina/Desktop/Otx2_3/Nuclear/KO_2i'
 %  '/Users/draina/Desktop/Otx2_3/Nuclear/WT_2i'
    '/Users/draina/Desktop/Otx2_3/Nuclear/KO_ESL'
 % '/Users/draina/Desktop/Otx2_3/Nuclear/WT_ESL'
    '/Users/draina/Desktop/Otx2_3/Nuclear/KO_N2C_12'
 %'/Users/draina/Desktop/Otx2_3/Nuclear/WT_N2C_12'
    '/Users/draina/Desktop/Otx2_3/Nuclear/KO_N2C_24'
 %'/Users/draina/Desktop/Otx2_3/Nuclear/WT_N2C_24'
    
    
    
    };

pathlist_labels = {
    %     'PDO3 (min ERK)'
    %     'PDO3 Release (max ERK)'
    %     'ES+Lif'
%     '1min'
%     '2min'
%     '3.5min'
%     '7.5min'
%     '15min'
%     '30min'
%     'ESL'
%     'PDO3'
    %'ESL'
    
    'KO_2i'
%    'WT_2i'
    'KO_ESL'
%    'WT_ESL'
    'KO_N2C_12'
%    'WT_N2C_12'
    'KO_N2C_24'
%    'WT_N2C_24'
    };


%Order Channel Labels according to image
ChannelLabel = {
    %         'DAPI'
    %         'NANOG'
    %         'SPRY4'
    %         'OCT3/4'
    
    'DAPI'
    'Otx2'
    'Trans'
    'Nanog'
    %
    %     'DAPI'
    %     'Trans'
    %     'Sensor'
    %     'ppErk'
    %
    
    };

ChannelCalcs = {                        %1 - Nuc, 2-Cyt, 3- Whole Cell, 4 - Nuc/Cyt ratio
    [0,0,0,0]
    [1,0,0,0]
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

boxplot = 0                                                            
conscat = 1
chanscat = 1


% Reorder Channels according to how you want to plot them (expects 4
% channels, so just duplicate or set one to garbage channel). Also expects
% Ch1 to be DAPI - doesn't quantify this channel!
ReorderChan = [ 1; 2; 4; 3];
%ReorderChan = [ 3; 2; 1; 4 ];
%ReorderChan = [ 1; 3; 4; 2 ];

inputTableType = 'csv';
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


for ctr2 = 1:length(pathlist)
    file(ctr2).path    = pathlist{ctr2};
    file(ctr2).name    = file(ctr2).path(find(file(ctr2).path==slashtype, 1, 'last')+1:end);
    file(ctr2).maindir = file(ctr2).path(1:(find(file(ctr2).path==slashtype, 1, 'last')-1));
    
    %Look for *.txt and .csv files
    dir_cell = [dir(fullfile(file(ctr2).path, '*.txt')); dir(fullfile(file(ctr2).path, '*.csv'))];
    
    %readtable is faster than loaddata!
    cnt1 = 1;
    cnt2 = 1;
    
    for c1 = 1:length(dir_cell)
        nuc_cyt_flag   = dir_cell(c1).name(1:3);
        
        if dir_cell(c1).name(end-2:end) == 'txt'
            delimiter      = '\t';
            temptable      = readtable([file(ctr2).path slashtype dir_cell(c1).name], 'delimiter', delimiter);
        elseif dir_cell(c1).name(end-2:end) == 'csv'
            delimiter      = 'comma';
            temptable      = readtable([file(ctr2).path slashtype dir_cell(c1).name], 'delimiter', delimiter);
        end
        
        switch nuc_cyt_flag
            case 'nuc'
                
                %Store Integrated Density
                intDenNuc_ch2{cnt1} = temptable.IntDen(temptable.Ch==ch2);
                intDenNuc_ch3{cnt1} = temptable.IntDen(temptable.Ch==ch3);
                intDenNuc_ch4{cnt1} = temptable.IntDen(temptable.Ch==ch4);
                
                %Store Area
                areaNuc_ch2{cnt1} = temptable.Area(temptable.Ch==ch2);
                areaNuc_ch3{cnt1} = temptable.Area(temptable.Ch==ch3);
                areaNuc_ch4{cnt1} = temptable.Area(temptable.Ch==ch4);
                
                %empty variables:
                intDenWhole_ch2{cnt1} = zeros(size(areaNuc_ch2{cnt1}));
                intDenWhole_ch3{cnt1} = zeros(size(areaNuc_ch2{cnt1}));
                intDenWhole_ch4{cnt1} = zeros(size(areaNuc_ch2{cnt1}));
                areaWhole_ch2{cnt1} = zeros(size(areaNuc_ch2{cnt1}));
                areaWhole_ch3{cnt1}   = zeros(size(areaNuc_ch2{cnt1}));
                areaWhole_ch4{cnt1}   = zeros(size(areaNuc_ch2{cnt1}));
                
                cnt1 = cnt1+1;
            case 'Cyt'
                
                %Store Integrated Density:
                intDenWhole_ch2{cnt2} = temptable.IntDen(temptable.Ch==ch2);
                intDenWhole_ch3{cnt2} = temptable.IntDen(temptable.Ch==ch3);
                intDenWhole_ch4{cnt2} = temptable.IntDen(temptable.Ch==ch4);
                
                %Store Area
                areaWhole_ch2{cnt2} = temptable.Area(temptable.Ch==ch2);
                areaWhole_ch3{cnt2} = temptable.Area(temptable.Ch==ch3);
                areaWhole_ch4{cnt2} = temptable.Area(temptable.Ch==ch4);
                
                %Empty Variables:
                intDenNuc_ch2{cnt2} = []; 
                intDenNuc_ch3{cnt2} = [];
                intDenNuc_ch4{cnt2} = [];
                areaNuc_ch2{cnt2}   = [];
                areaNuc_ch3{cnt2}   = [];
                areaNuc_ch4{cnt2}   = [];
                
                cnt2 = cnt2+1;
        end
    end
    clear cnt1 cnt2
    
    
    %Concat. everything from loop
    totIntWhole_2 = cell2mat(intDenWhole_ch2');
    totIntWhole_3 = cell2mat(intDenWhole_ch3');
    totIntWhole_4 = cell2mat(intDenWhole_ch4');
    
    totAreaWhole_2 = cell2mat(areaWhole_ch2');
    totAreaWhole_3 = cell2mat(areaWhole_ch3');
    totAreaWhole_4 = cell2mat(areaWhole_ch4');
    
    totIntNuc_2 = cell2mat(intDenNuc_ch2');
    totIntNuc_3 = cell2mat(intDenNuc_ch3');
    totIntNuc_4 = cell2mat(intDenNuc_ch4');
    
    totAreaNuc_2 = cell2mat(areaNuc_ch2');
    totAreaNuc_3 = cell2mat(areaNuc_ch3');
    totAreaNuc_4 = cell2mat(areaNuc_ch4');
    
    
    %Estimating Cytoplasm:
    
    totIntCyt_2 = totIntWhole_2-totIntNuc_2;
    totIntCyt_3 = totIntWhole_3-totIntNuc_3;
    totIntCyt_4 = totIntWhole_4-totIntNuc_4;
    
    totAreaCyt_2 = totAreaWhole_2-totAreaNuc_2;
    totAreaCyt_3 = totAreaWhole_3-totAreaNuc_3;
    totAreaCyt_4 = totAreaWhole_4-totAreaNuc_4;
    
    
    %Calculating Means:
    meanNuc_2 = totIntNuc_2./totAreaNuc_2;
    meanNuc_3 = totIntNuc_3./totAreaNuc_3;
    meanNuc_4 = totIntNuc_4./totAreaNuc_4;
    
    meanCyt_2 = totIntCyt_2./totAreaCyt_2;
    meanCyt_3 = totIntCyt_3./totAreaCyt_3;
    meanCyt_4 = totIntCyt_4./totAreaCyt_4;
    
    meanWhole_2 = totIntWhole_2./totAreaWhole_2;
    meanWhole_3 = totIntWhole_3./totAreaWhole_3;
    meanWhole_4 = totIntWhole_4./totAreaWhole_4;
    
    
    %Calculate N/C Ratio:
    meanRatio_2 = meanNuc_2./meanCyt_2;
    meanRatio_3 = meanNuc_3./meanCyt_3;
    meanRatio_4 = meanNuc_4./meanCyt_4;
    
    
    %Assign the Channel Indices for error-free plot labeling:
    meanNuc_2(1:length(meanNuc_2),2)=2;
    meanNuc_3(1:length(meanNuc_3),2)=3;
    meanNuc_4(1:length(meanNuc_4),2)=4;
    
    meanCyt_2(1:length(meanCyt_2),2)=2;
    meanCyt_3(1:length(meanCyt_3),2)=3;
    meanCyt_4(1:length(meanCyt_4),2)=4;
    
    meanWhole_2(1:length(meanWhole_2),2)=2;
    meanWhole_3(1:length(meanWhole_3),2)=3;
    meanWhole_4(1:length(meanWhole_4),2)=4;
    
    meanRatio_2(1:length(meanRatio_2),2)=2;
    meanRatio_3(1:length(meanRatio_3),2)=3;
    meanRatio_4(1:length(meanRatio_4),2)=4;
    
    
    %%
    
    %%%%%%%%%%%%%%%%%
    %Saving desired data
    resvec_calc2{1,ctr2} = meanNuc_2;   %For Channel2
    resvec_calc2{2,ctr2} = meanCyt_2;
    resvec_calc2{3,ctr2} = meanWhole_2;
    resvec_calc2{4,ctr2} = meanRatio_2;
    
    resvec_calc3{1,ctr2} = meanNuc_3;   %For Channel3
    resvec_calc3{2,ctr2} = meanCyt_3;
    resvec_calc3{3,ctr2} = meanWhole_3;
    resvec_calc3{4,ctr2} = meanRatio_3;
    
    resvec_calc4{1,ctr2} = meanNuc_4;      %For Channel4
    resvec_calc4{2,ctr2} = meanCyt_4;
    resvec_calc4{3,ctr2} = meanWhole_4;
    resvec_calc4{4,ctr2} = meanRatio_4;
    
    resvec_another{1, ctr2} = [meanNuc_2; meanNuc_3; meanNuc_4];
    resvec_another{2, ctr2} = [meanCyt_2; meanCyt_3; meanCyt_4];
    resvec_another{3, ctr2} = [meanWhole_2; meanWhole_3; meanWhole_4];
    resvec_another{4, ctr2} = [meanRatio_2; meanRatio_3; meanRatio_4];
    %%%%%%%%%%%%%%%%%
    
    
    %Clear Vars for loop:
    clear meanNuc_2 meanNuc_3 meanNuc_4 meanCyt_2 meanCyt_3 meanCyt_4 meanWhole_2 meanWhole_3 meanWhole_4 meanRatio_2 meanRatio_3 meanRatio_4 temptable ...
        intDenNuc_ch2 intDenNuc_ch3 intDenNuc_ch4 areaNuc_ch2 areaNuc_ch3 areaNuc_ch4 ...
        intDenWhole_ch2 intDenWhole_ch3 intDenWhole_ch4 areaWhole_ch2 areaWhole_ch3 areaWhole_ch4
end


if boxplot==1
    
    %% Per-Channel Box Plots (i.e. across all treatment groups for one channel)
    activeChans = cellfun(@(x) sum(x),ChannelCalcsR)>0; %Reads the ChannelCalcs variable to decide which channels need calcs
    plotflag = [1 0 0];
    scatx = 0;
    scaty = 0;
    for rr = 1:length(activeChans) %Channels
        if activeChans(rr)>0
            for cc = 1:length(ChannelCalcsR{rr}) %Calculations
                calcType = ChannelCalcsR{rr}(cc);
                if calcType>0
                    
                    switch rr
                        case 2
                            resvec = resvec_calc2(calcType,:);
                        case 3
                            resvec = resvec_calc3(calcType,:);
                        case 4
                            resvec = resvec_calc4(calcType,:);
                    end
                    
                    chlabel = ChannelLabel(ReorderChan(rr));
                    calclbl2 = char(calclbl{calcType});
                    IF_ncplot(plotflag, resvec,scatx, scaty, pathlist_labels, chlabel, calclbl2, file, slashtype)
                end
            end
        end
    end
end

if chanscat ==1
    
    %% Per Channel Scatter
    clear chlabel
    xchan = 2;
    ychan = 3;
    %calcType : 1-Nuc, 2-Cyt, 3-WholeCell, 4-N/C Ratio
    xcalcType = 1;
    ycalcType = 1;
    
    plotflag = [0 1 0];
    chlabel{1,1} = char(ChannelLabel(ReorderChan(xchan)));
    chlabel{1,2} = char(ChannelLabel(ReorderChan(ychan)));
    
    chlabel{2,1} = char(calclbl{xcalcType});
    chlabel{2,2} = char(calclbl{ycalcType});
    resvec = 0;
    for nn = 1:length(file)
        switch xchan
            case(2)
                scatx = resvec_calc2{xcalcType,nn}(:,1);
            case(3)
                scatx = resvec_calc3{xcalcType,nn}(:,1);
            case(4)
                scatx = resvec_calc4{xcalcType,nn}(:,1);
        end
        switch ychan
            case(2)
                scaty = resvec_calc2{ycalcType,nn}(:,1);
            case(3)
                scaty = resvec_calc3{ycalcType,nn}(:,1);
            case(4)
                scaty = resvec_calc4{ycalcType,nn}(:,1);
        end

        plotflag(2,1) = nn;
        IF_ncplot(plotflag, resvec,scatx, scaty, char(pathlist_labels{nn}), chlabel, calclbl, file, slashtype)
    end
    
end


if conscat ==1
    %% Consolidated Scatter:
    clear chlabel
    TreatmentList = [1 2];
    xchan = 2;
    ychan = 3;
    %calcType : 1-Nuc, 2-Cyt, 3-WholeCell, 4-N/C Ratio
    xcalcType = 1;
    ycalcType = 1;
    
    plotflag = [0 0 1];
    chlabel{1} = [char(ChannelLabel(ReorderChan(xchan))) ' ' char(calclbl(xcalcType))];
    chlabel{2} = [char(ChannelLabel(ReorderChan(ychan))) ' ' char(calclbl(ycalcType))];
    resvec = 0;
    
    switch xchan
        case(2)
            scatx = resvec_calc2(xcalcType,:);
        case(3)
            scatx = resvec_calc3(xcalcType,:);
        case(4)
            scatx = resvec_calc4(xcalcType,:);
    end
    switch ychan
        case(2)
            scaty = resvec_calc2(ycalcType,:);
        case(3)
            scaty = resvec_calc3(ycalcType,:);
        case(4)
            scaty = resvec_calc4(ycalcType,:);
    end
    plotflag(3,1:length(TreatmentList)) = TreatmentList;
    plotflag(2,1) = length(file); %For colours
    IF_ncplot(plotflag, resvec,scatx, scaty, pathlist_labels, chlabel, calclbl, file, slashtype)
    
    
    
end
