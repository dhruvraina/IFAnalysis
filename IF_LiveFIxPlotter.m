%IF_Live plotter
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



boxplot  = 1;
chanscat = 1;
conscat  = 1;
plotflag.singScat = 0;
plotflag.conScat  = 0;
plotflag.boxplot  = 0;

%% User Inputs:
pathlist = {
    
'/Users/draina/Desktop/2018_E1_3_otx2Fixed/Nuclear Intensities/Only Tracked Cells/KO_005'
'/Users/draina/Desktop/2018_E1_3_otx2Fixed/Nuclear Intensities/Only Tracked Cells/KO_100'
};

pathlist_labels = {
    'KO 005ng'
    'KO 100ng'
    };


%% Main

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
            keyboard
            %Put whatever else here, just change the 'delimiter'
        elseif dir_cell(c1).name(end-2:end) == 'csv'
            delimiter      = ';';
            temptable      = readtable([file(ctr2).path slashtype dir_cell(c1).name], 'delimiter', delimiter);
        end
        
        mean_ch1{cnt1}      = temptable.MeanNanog;
        mean_ch2{cnt1}      = temptable.MeanOtx2;
        mean_actfrac{cnt1}  = temptable.ActiveFraction;
        mean_peakfreq{cnt1} = temptable.PeakFrequency;
    end
    clear cnt1 cnt2
    
    totInt_1 = cell2mat(mean_ch1');
    totInt_2 = cell2mat(mean_ch2');
    totFrac = cell2mat(mean_actfrac');
    totFreq = cell2mat(mean_peakfreq'); 
    
    %Calculate N/C Ratio:
    meanRatio = totInt_1./totInt_2;
        
    %Assign the Channel Indices for error-free plot labeling:
    totInt_1(1:length(totInt_1),2)=1;
    totInt_2(1:length(totInt_2),2)=2;
    totFrac(1:length(totFrac),2)=ctr2;
    totFreq(1:length(totFreq),2)=ctr2;
    meanRatio(1:length(meanRatio),2)=ctr2;
    
    %%
    
    %%%%%%%%%%%%%%%%%
    %Saving desired data
    resvec_calc{1,ctr2} = totInt_1;   
    resvec_calc{2,ctr2} = totInt_2;
    resvec_calc{3,ctr2} = totFrac;
    resvec_calc{4,ctr2} = totFreq;
    resvec_calc{5,ctr2} = meanRatio;

    %%%%%%%%%%%%%%%%%
    
    
    %Clear Vars for loop:
    clear temptable totInt_1 totInt_2 totFrac totFreq meanRatio
end

    calcLabels= {'MeanNanog'...
                 'MeanOtx2'...
                 'ActiveFraction'...
                 'PeakFrequency'...
                 'NanogdivOtx2'};
             

if boxplot==1
    
    %% Box Plots
    plotflag.type = 'boxplot';
    scatx = 0;
    scaty = 0;
    calcType = [1,2,3,4,5];
            for cc = calcType %Calculations
                if cc>0
                    chlabel = 'null';
                    resvec = resvec_calc(cc,:);
                    calclbl2 = char(calcLabels(cc));
                    IF_ncplot(plotflag, resvec,scatx, scaty, pathlist_labels, chlabel, calclbl2, file, slashtype)
                end
            end
    
end

    calcLabels= {'MeanNanog'...
                 'MeanOtx2'...
                 'ActiveFraction'...
                 'PeakFrequency'...
                 'NanogdivOtx2'};

if chanscat ==1
    
    %% Per Channel Scatter
    clear chlabel
    calclbl = ''; 
    xchan = 5;
    ychan = 4;
    %calcType : 1-Nuc, 2-Cyt, 3-WholeCell, 4-N/C Ratio
    xcalcType = '_';
    ycalcType = '_';
    
    plotflag.type = 'SingleScatter';
    plotflag.corrprint =1;
    chlabel{1,1} = char(calcLabels(xchan));
    chlabel{1,2} = char(calcLabels(ychan));
    
    chlabel{2,1} = char(xcalcType);
    chlabel{2,2} = char(ycalcType);
    
    resvec = 0;
    for nn = 1:length(file)

                scatx = resvec_calc{xchan,nn}(:,1);
                scaty = resvec_calc{ychan,nn}(:,1);

        plotflag.singTreat = nn;
        IF_ncplot(plotflag, resvec,scatx, scaty, char(pathlist_labels{nn}), chlabel, calclbl, file, slashtype)
    end
    
end


if conscat ==1
    %% Consolidated Scatter:
    clear chlabel
    TreatmentList = [1 2];
    xchan = 5;
    ychan = 3;
    %calcType : 1-Nuc, 2-Cyt, 3-WholeCell, 4-N/C Ratio
    xcalcType = '_';
    ycalcType = '_';
    
    plotflag.type = 'ConScatter'; 
    chlabel{1} = char(calcLabels(xchan));
    chlabel{2} = char(calcLabels(ychan));
    
    resvec = 0;
    
    scatx = resvec_calc(xchan,:);
    scaty = resvec_calc(ychan,:);

    plotflag.conTreat = TreatmentList;
    plotflag.conColour = length(file); %For colours
    IF_ncplot(plotflag, resvec,scatx, scaty, pathlist_labels, chlabel, calclbl, file, slashtype)
    
    
    
end
