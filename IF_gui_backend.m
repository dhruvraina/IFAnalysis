%IF Plotter
%Description: Backend for IF_GUI. Plots .csv or .txt data by grouping
%a table of values by the header 'Ch' into channels.
%Last Edit: 180822
%Author: draina


function IF_gui_backend(file, inputs, calcs, outputs)
%% Main
%Can be Reordered here:
ch1 = 1;
ch2 = 2;
ch3 = 3;
ch4 = 4;
plotflag.imageFormat = outputs.imageFormat;
plotflag.margDist    = outputs.margDist;

for ctr2 = 1:length(file.treatmentfold)
    
    file.tdir = [file.path file.slashtype file.treatmentfold{ctr2}];
    
    %Look for *.txt and .csv files
    dir_cell = [dir(fullfile(file.tdir, '*.txt')); dir(fullfile(file.tdir, '*.csv'))];
    
    %Error Handling:
    if isempty(dir_cell)
        errordlg('Check directory names, error in file struct')
    end
    
    cnt1 = 1;
    cnt2 = 1;
    
    %% 2. Concatenate all data from individual .txt or .csv
    for c1 = 1:length(dir_cell)
        
        
        %Automatically detect if .csv or .txt and apply readtable
        if dir_cell(c1).name(end-2:end)     == 'txt'
            delimiter      = '\t';
            temptable      = readtable([file.tdir file.slashtype dir_cell(c1).name], 'delimiter', delimiter);
        elseif dir_cell(c1).name(end-2:end) == 'csv'
            delimiter      = 'comma';
            temptable      = readtable([file.tdir file.slashtype dir_cell(c1).name], 'delimiter', delimiter);
        end
        
        
        %Read in prefix to determine if 'Nuc' or 'Cyt'
        nuc_cyt_flag   = dir_cell(c1).name(1:3);
        
        
        %Error Handling:
        if size(temptable,2)<11
            errordlg('Number of columns in .txt file is incorrect. Please use *RAW* data from imageJ script')
        end
        
        switch nuc_cyt_flag
            case inputs.nucMaskPrefix
                %%%%%%%%%%%%----------------------------------------------_%%%%%%-_%_%_%
                %Store Integrated Density
                intDenNuc_ch1{cnt1} = temptable.IntDen(temptable.Ch==ch1);
                intDenNuc_ch2{cnt1} = temptable.IntDen(temptable.Ch==ch2);
                intDenNuc_ch3{cnt1} = temptable.IntDen(temptable.Ch==ch3);
                intDenNuc_ch4{cnt1} = temptable.IntDen(temptable.Ch==ch4);
                
                %Store Area
                areaNuc_ch1{cnt1} = temptable.Area(temptable.Ch==ch1);
                areaNuc_ch2{cnt1} = temptable.Area(temptable.Ch==ch2);
                areaNuc_ch3{cnt1} = temptable.Area(temptable.Ch==ch3);
                areaNuc_ch4{cnt1} = temptable.Area(temptable.Ch==ch4);
                
                %empty variables:
                [intDenWhole_ch1{cnt1}, intDenWhole_ch2{cnt1}, ...
                    intDenWhole_ch3{cnt1}, intDenWhole_ch4{cnt1}, ...
                    areaWhole_ch1{cnt1},   areaWhole_ch2{cnt1}, ...
                    areaWhole_ch3{cnt1},   areaWhole_ch4{cnt1}, ...
                    ] = deal(zeros(size(areaNuc_ch2{cnt1})));
                
                cnt1 = cnt1+1;
            case inputs.cytMaskPrefix
                
                %Store Integrated Density:
                intDenWhole_ch1{cnt2} = temptable.IntDen(temptable.Ch==ch1);
                intDenWhole_ch2{cnt2} = temptable.IntDen(temptable.Ch==ch2);
                intDenWhole_ch3{cnt2} = temptable.IntDen(temptable.Ch==ch3);
                intDenWhole_ch4{cnt2} = temptable.IntDen(temptable.Ch==ch4);
                
                %Store Area
                areaWhole_ch1{cnt2} = temptable.Area(temptable.Ch==ch1);
                areaWhole_ch2{cnt2} = temptable.Area(temptable.Ch==ch2);
                areaWhole_ch3{cnt2} = temptable.Area(temptable.Ch==ch3);
                areaWhole_ch4{cnt2} = temptable.Area(temptable.Ch==ch4);
                
                %Empty Variables:
                [intDenNuc_ch1{cnt2}, intDenNuc_ch2{cnt2}, ...
                    intDenNuc_ch3{cnt2}, intDenNuc_ch4{cnt2}, ...
                    areaNuc_ch1{cnt2},   areaNuc_ch2{cnt2}, ...
                    areaNuc_ch3{cnt2},   areaNuc_ch4{cnt2}, ...
                    ] = deal([]);
                
                cnt2 = cnt2+1;
        end
    end
    clear cnt1 cnt2
    
    
    %Concat. everything from loop
    totIntWhole_1 = cell2mat(intDenWhole_ch1');
    totIntWhole_2 = cell2mat(intDenWhole_ch2');
    totIntWhole_3 = cell2mat(intDenWhole_ch3');
    totIntWhole_4 = cell2mat(intDenWhole_ch4');
    
    totAreaWhole_1 = cell2mat(areaWhole_ch1');
    totAreaWhole_2 = cell2mat(areaWhole_ch2');
    totAreaWhole_3 = cell2mat(areaWhole_ch3');
    totAreaWhole_4 = cell2mat(areaWhole_ch4');
    
    totIntNuc_1 = cell2mat(intDenNuc_ch1');
    totIntNuc_2 = cell2mat(intDenNuc_ch2');
    totIntNuc_3 = cell2mat(intDenNuc_ch3');
    totIntNuc_4 = cell2mat(intDenNuc_ch4');
    
    totAreaNuc_1 = cell2mat(areaNuc_ch1');
    totAreaNuc_2 = cell2mat(areaNuc_ch2');
    totAreaNuc_3 = cell2mat(areaNuc_ch3');
    totAreaNuc_4 = cell2mat(areaNuc_ch4');
    
    
    %Estimating Cytoplasm:
    totIntCyt_1 = totIntWhole_1-totIntNuc_1;
    totIntCyt_2 = totIntWhole_2-totIntNuc_2;
    totIntCyt_3 = totIntWhole_3-totIntNuc_3;
    totIntCyt_4 = totIntWhole_4-totIntNuc_4;
    
    totAreaCyt_1 = totAreaWhole_1-totAreaNuc_1;
    totAreaCyt_2 = totAreaWhole_2-totAreaNuc_2;
    totAreaCyt_3 = totAreaWhole_3-totAreaNuc_3;
    totAreaCyt_4 = totAreaWhole_4-totAreaNuc_4;
    
    
    %Calculating Means:
    meanNuc_1 = totIntNuc_1./totAreaNuc_1;
    meanNuc_2 = totIntNuc_2./totAreaNuc_2;
    meanNuc_3 = totIntNuc_3./totAreaNuc_3;
    meanNuc_4 = totIntNuc_4./totAreaNuc_4;
    
    meanCyt_1 = totIntCyt_1./totAreaCyt_1;
    meanCyt_2 = totIntCyt_2./totAreaCyt_2;
    meanCyt_3 = totIntCyt_3./totAreaCyt_3;
    meanCyt_4 = totIntCyt_4./totAreaCyt_4;
    
    meanWhole_1 = totIntWhole_1./totAreaWhole_1;
    meanWhole_2 = totIntWhole_2./totAreaWhole_2;
    meanWhole_3 = totIntWhole_3./totAreaWhole_3;
    meanWhole_4 = totIntWhole_4./totAreaWhole_4;
    
    
    %Calculate N/C Ratio:
    meanRatio_1 = meanNuc_1./meanCyt_1;
    meanRatio_2 = meanNuc_2./meanCyt_2;
    meanRatio_3 = meanNuc_3./meanCyt_3;
    meanRatio_4 = meanNuc_4./meanCyt_4;
    
    
    %Assign the Channel Indices for error-free plot labeling:
    meanNuc_1(1:length(meanNuc_1),2)=1;
    meanNuc_2(1:length(meanNuc_2),2)=2;
    meanNuc_3(1:length(meanNuc_3),2)=3;
    meanNuc_4(1:length(meanNuc_4),2)=4;
    
    meanCyt_1(1:length(meanCyt_1),2)=1;
    meanCyt_2(1:length(meanCyt_2),2)=2;
    meanCyt_3(1:length(meanCyt_3),2)=3;
    meanCyt_4(1:length(meanCyt_4),2)=4;
    
    meanWhole_1(1:length(meanWhole_1),2)=1;
    meanWhole_2(1:length(meanWhole_2),2)=2;
    meanWhole_3(1:length(meanWhole_3),2)=3;
    meanWhole_4(1:length(meanWhole_4),2)=4;
    
    meanRatio_1(1:length(meanRatio_1),2)=1;
    meanRatio_2(1:length(meanRatio_2),2)=2;
    meanRatio_3(1:length(meanRatio_3),2)=3;
    meanRatio_4(1:length(meanRatio_4),2)=4;
    
    
    %%
    
    %%%%%%%%%%%%%%%%%
    %Saving desired data
    %Channel1
    resvec_calc1{1,ctr2} = meanNuc_1;
    resvec_calc1{2,ctr2} = meanCyt_1;
    resvec_calc1{3,ctr2} = meanWhole_1;
    resvec_calc1{4,ctr2} = meanRatio_1;
    resvec_calc1{5,ctr2} = meanRatio_1./1;
    
    %Channel2
    resvec_calc2{1,ctr2} = meanNuc_2;
    resvec_calc2{2,ctr2} = meanCyt_2;
    resvec_calc2{3,ctr2} = meanWhole_2;
    resvec_calc2{4,ctr2} = meanRatio_2;
    resvec_calc2{5,ctr2} = meanRatio_2./1;
    
    %Channel3
    resvec_calc3{1,ctr2} = meanNuc_3;
    resvec_calc3{2,ctr2} = meanCyt_3;
    resvec_calc3{3,ctr2} = meanWhole_3;
    resvec_calc3{4,ctr2} = meanRatio_3;
    resvec_calc3{5,ctr2} = meanRatio_3./1;
    
    %Channel4
    resvec_calc4{1,ctr2} = meanNuc_4;
    resvec_calc4{2,ctr2} = meanCyt_4;
    resvec_calc4{3,ctr2} = meanWhole_4;
    resvec_calc4{4,ctr2} = meanRatio_4;
    resvec_calc4{5,ctr2} = meanRatio_4./1;
    
    %BackupCalcs
    resvec_another{1, ctr2} = [meanNuc_1; meanNuc_2; meanNuc_3; meanNuc_4];
    resvec_another{2, ctr2} = [meanCyt_1; meanCyt_2; meanCyt_3; meanCyt_4];
    resvec_another{3, ctr2} = [meanWhole_1; meanWhole_2; meanWhole_3; meanWhole_4];
    resvec_another{4, ctr2} = [meanRatio_1; meanRatio_2; meanRatio_3; meanRatio_4];
    resvec_another{5, ctr2} = [meanRatio_1./1; meanRatio_2./1; meanRatio_3./1; meanRatio_4./1];
    %%%%%%%%%%%%%%%%%
    
    
    %Clear Vars for loop:
    clear meanNuc_1      meanNuc_2       meanNuc_3       meanNuc_4 ...
        meanCyt_1      meanCyt_2       meanCyt_3       meanCyt_4 ...
        meanWhole_1    meanWhole_2     meanWhole_3     meanWhole_4 ...
        meanRatio_1    meanRatio_2     meanRatio_3     meanRatio_4 temptable ...
        intDenNuc_ch1  intDenNuc_ch2   intDenNuc_ch3   intDenNuc_ch4 ...
        areaNuc_ch1    areaNuc_ch2     areaNuc_ch3     areaNuc_ch4 ...
        intDenWhole_ch1 intDenWhole_ch2 intDenWhole_ch3 intDenWhole_ch4 ...
        areaWhole_ch1  areaWhole_ch2   areaWhole_ch3    areaWhole_ch4
end

%% |----------- NORMALIZATION ------------|
if calcs.normFlag
    switch calcs.normType
        case 'max'
            normvals = cellfun(@(x) max(x(:,1)), resvec_calc1(:,calcs.normTo));
        case 'min'
            normvals = cellfun(@(x) min(x(:,1)), resvec_calc1(:,calcs.normTo));
    end
    
    norm_resvec_calc1 = [];
    norm_resvec_calc2 = [];
    norm_resvec_calc3 = [];
    norm_resvec_calc4 = [];
    
    for c2 = 1:size(resvec_calc1,1)
        norm_resvec_calc1 = [norm_resvec_calc1; cellfun(@(x) x./normvals(c2), resvec_calc1(c2,:), 'UniformOutput', 0)];
        norm_resvec_calc2 = [norm_resvec_calc2; cellfun(@(x) x./normvals(c2), resvec_calc2(c2,:), 'UniformOutput', 0)];
        norm_resvec_calc3 = [norm_resvec_calc3; cellfun(@(x) x./normvals(c2), resvec_calc3(c2,:), 'UniformOutput', 0)];
        norm_resvec_calc4 = [norm_resvec_calc4; cellfun(@(x) x./normvals(c2), resvec_calc4(c2,:), 'UniformOutput', 0)];
    end
end


%% --------------  A. Boxplots  --------:

if outputs.boxplot==1
    
    %Swap save variables during normalization:
    if calcs.normFlag
        resvec_ch1 = norm_resvec_calc1;
        resvec_ch2 = norm_resvec_calc2;
        resvec_ch3 = norm_resvec_calc3;
        resvec_ch4 = norm_resvec_calc4;
        
    else
        resvec_ch1 = resvec_calc1;
        resvec_ch2 = resvec_calc2;
        resvec_ch3 = resvec_calc3;
        resvec_ch4 = resvec_calc4;
        
    end
    
    
    activeChans = find(inputs.activeChans);                                %Reads the ChannelCalcs variable to decide which channels need calcs
    if ~isempty(activeChans)
        
        %Inputs for IF_ncplot.m
        plotflag.type = 'boxplot2'                                         %boxplot1 is the basic NotBoxPlot;
        scatx = 0;                                                         %boxplot2 is the UnivarScatter plot;
        scaty = 0;
        
        for cc = activeChans                                               %Channels - in columns
            for rr = 1:size(calcs.all,1)                                   %Loop through different calculations - in rows
                
                if calcs.all(rr,cc)>0
                    
                    switch cc
                        case 1
                            resvec = resvec_ch1(rr,:);
                        case 2
                            resvec = resvec_ch2(rr,:);
                        case 3
                            resvec = resvec_ch3(rr,:);
                        case 4
                            resvec = resvec_ch4(rr,:);
                    end
                    
                    
                    switch outputs.plotMode
                        case 'subset'
                            resvec = resvec(:,outputs.limconscatID);
                            treatmentLabels = file.treatmentLabels(outputs.limconscatID);
                        case 'all'
                            treatmentLabels = file.treatmentLabels;
                    end
                    
                    %Calculate Limits
                    if outputs.boxAutoY
                        lims.boxmax =  ceil(max(cellfun(@(x) max(x(:,1)), resvec))/10)*10;
                        lims.boxmin = floor(min(cellfun(@(x) min(x(:,1)), resvec))/10)*10;
                    else
                        lims.boxmin = outputs.boxYlim(1);
                        lims.boxmax = outputs.boxYlim(2);
                    end
                    
                    lims.boxmean= cellfun(@(x) mean(x(:,1)), resvec); %mainly for printing
                    chlabel  = inputs.ChannelLabel{cc};
                    calclbl2 = calcs.label{cc};
                    IF_ncplot(plotflag, resvec,scatx, scaty, treatmentLabels, chlabel, calclbl2, file, lims)
                end
            end
        end
        
        
    end
end


%% --------------  B. Single Treatment Scatter Plots  --------:
if outputs.chanscat ==1
    
    %Swap save variables during normalization:
    if calcs.normFlag
        resvec_ch1 = norm_resvec_calc1;
        resvec_ch2 = norm_resvec_calc2;
        resvec_ch3 = norm_resvec_calc3;
        resvec_ch4 = norm_resvec_calc4;
        
    else
        resvec_ch1 = resvec_calc1;
        resvec_ch2 = resvec_calc2;
        resvec_ch3 = resvec_calc3;
        resvec_ch4 = resvec_calc4;
        
    end
    
    
    clear chlabel
    xchan  = outputs.scatX;
    ychan  = outputs.scatY;
    lims.x = outputs.scatXlim;
    lims.y = outputs.scatYlim;
    
    %calcType : 1-Nuc, 2-Cyt, 3-WholeCell, 4-N/C Ratio, 5-C/N Ratio
    xcalcType = find(calcs.all(:,xchan));
    ycalcType = find(calcs.all(:,ychan));
    
    plotflag.type      = 'SingleScatter';
    plotflag.corrprint = outputs.CorrLine;
    
    chlabel{1,1} = inputs.ChannelLabel{xchan};
    chlabel{1,2} = inputs.ChannelLabel{ychan};
    chlabel{2,1} = calcs.label{xchan};
    chlabel{2,2} = calcs.label{ychan};
    
    %Unnecessary Inputs:
    resvec  = 0;
    calclbl = calcs.label;
    
    %Subset of treatments to plot:
    switch outputs.plotMode
        case 'subset'
            treatmentList = outputs.limconscatID;
        case 'all'
            treatmentList = 1:length(file.treatmentfold);
    end
    
    
    %Assigning X and Y vectors
    for nn = treatmentList
        switch xchan
            case(1)
                scatx    = resvec_ch1{xcalcType,nn}(:,1);
                resvecX   = resvec_ch1;
            case(2)
                scatx = resvec_ch2{xcalcType,nn}(:,1);
                resvecX   = resvec_ch2;
            case(3)
                scatx = resvec_ch3{xcalcType,nn}(:,1);
                resvecX   = resvec_ch3;
            case(4)
                scatx = resvec_ch4{xcalcType,nn}(:,1);
                resvecX   = resvec_ch4;
        end
        
        switch ychan
            case(1)
                scaty = resvec_ch1{ycalcType,nn}(:,1);
                resvecY   = resvec_ch1;
            case(2)
                scaty = resvec_ch2{ycalcType,nn}(:,1);
                resvecY   = resvec_ch2;
            case(3)
                scaty = resvec_ch3{ycalcType,nn}(:,1);
                resvecY   = resvec_ch3;
            case(4)
                scaty = resvec_ch4{ycalcType,nn}(:,1);
                resvecY   = resvec_ch4;
        end
        
        
        %Set auto limits:
        if outputs.scatAutoX
            lims.x = [floor(min(cellfun(@(x) min(x(:,1)), resvecX(xcalcType,:)))/10)*10 ...
                ceil(max(cellfun(@(x) max(x(:,1)), resvecX(xcalcType,:)))/10)*10];
        end
        
        if outputs.scatAutoY
            lims.y = [floor(min(cellfun(@(x) min(x(:,1)), resvecY(xcalcType,:)))/10)*10 ...
                ceil(max(cellfun(@(x) max(x(:,1)), resvecY(xcalcType,:)))/10)*10];
        end
        
        %Plotflag keeps track of the colours for single scatters
        plotflag.singTreat = nn;
        IF_ncplot(plotflag, resvec, scatx, scaty, file.treatmentLabels{nn}, chlabel, calclbl, file, lims)
    end
    
end




%% --------------  C. Multiple Treatments In One Scatter  --------:
if outputs.conscat==1
    
    %Swap save variables during normalization:
    if calcs.normFlag
        resvec_ch1 = norm_resvec_calc1;
        resvec_ch2 = norm_resvec_calc2;
        resvec_ch3 = norm_resvec_calc3;
        resvec_ch4 = norm_resvec_calc4;
        
    else
        resvec_ch1 = resvec_calc1;
        resvec_ch2 = resvec_calc2;
        resvec_ch3 = resvec_calc3;
        resvec_ch4 = resvec_calc4;
        
    end
    
    
    clear chlabel
    plotflag.type     = 'ConScatter';
    
    switch outputs.plotMode
        case 'subset'
            plotflag.conTreat = outputs.limconscatID;
        case 'all'
            plotflag.conTreat = 1:length(file.treatmentLabels);                    %Replace this with the actual list of treatments
    end
    lims.x = outputs.scatXlim;
    lims.y = outputs.scatYlim;
    xchan  = outputs.scatX;
    ychan  = outputs.scatY;
    
    %calcType : 1-Nuc, 2-Cyt, 3-WholeCell, 4-N/C Ratio 5-C/N Ratio
    xcalcType = find(calcs.all(:,xchan));
    ycalcType = find(calcs.all(:,ychan));
    
    
    chlabel{1} = [inputs.ChannelLabel{xchan} ' ' calcs.label{xchan}];
    chlabel{2} = [inputs.ChannelLabel{ychan} ' ' calcs.label{ychan}];
    
    %Unnecessary Inputs:
    resvec  = 0;
    calclbl = calcs.label;
    
    %Reading out X and Y into vectors:
    switch xchan
        case(1)
            scatx = resvec_ch1(xcalcType,:);
            resvecX   = resvec_ch1;
        case(2)
            scatx = resvec_ch2(xcalcType,:);
            resvecX   = resvec_ch2;
        case(3)
            scatx = resvec_ch3(xcalcType,:);
            resvecX   = resvec_ch3;
        case(4)
            scatx = resvec_ch4(xcalcType,:);
            resvecX   = resvec_ch4;
    end
    
    switch ychan
        case(1)
            scaty = resvec_ch1(ycalcType,:);
            resvecY   = resvec_ch1;
        case(2)
            scaty = resvec_ch2(ycalcType,:);
            resvecY   = resvec_ch2;
        case(3)
            scaty = resvec_ch3(ycalcType,:);
            resvecY   = resvec_ch3;
        case(4)
            scaty = resvec_ch4(ycalcType,:);
            resvecY   = resvec_ch4;
    end
    
    %Keep colours consistent between ConScat and SingleScat:
    plotflag.colours = 1:length(file.treatmentLabels);
    
    %Auto limits
    if outputs.scatAutoX
        lims.x = [floor(min(cellfun(@(x) min(x(:,1)), resvecX(xcalcType,:)))/10)*10 ...
            ceil(max(cellfun(@(x) max(x(:,1)), resvecX(xcalcType,:)))/10)*10];
    end
    
    if outputs.scatAutoY
        lims.y = [floor(min(cellfun(@(x) min(x(:,1)), resvecY(xcalcType,:)))/10)*10 ...
            ceil(max(cellfun(@(x) max(x(:,1)), resvecY(xcalcType,:)))/10)*10];
    end
    IF_ncplot(plotflag, resvec,scatx, scaty, file.treatmentLabels, chlabel, calclbl, file, lims)
    
    
end
save([file.outpath file.slashtype char(cleanNames({file.experimentName}, '_')) 'results.mat'], 'resvec_calc1', 'resvec_calc2', 'resvec_calc3', 'resvec_calc4');


calcType = 1
r_name = [];
if ~isempty(r_name)
    fileID = fopen([file.outpath file.slashtype '_Reformat.txt'],'w');
    for tnum = 1:size(treatmentLabels,2)
        tcell = horzcat(resvec_calc1{calcType, tnum}(:,1), resvec_calc2{calcType, tnum}(:,1), resvec_calc3{calcType, tnum}(:,1), resvec_calc4{calcType, tnum}(:,1));
        fprintf(fileID, '%1$s\t\t\t\t\t\r\n', treatmentLabels{tnum})
        
        fprintf(fileID, '%1$s\t %2$s\t %3$s\t %4$s\t', calcs.label(1))
    end
    
    
    fprintf(fileID,'%1$s\t %2$s\t %3$s\t %4$s\t %5$s\t %6$s\t %7$s\t %8$s\t %9$s\r\n','TrackID', 'RiseTimes', 'FallTimes', 'ActiveFreq', 'ActiveLength', 'TotalLength', 'ActiveFraction', 'PeakNums', 'PeakFrequency');
    fprintf(fileID, '%d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\r\n', horzcat(peaknums{cc,ii}(:,2), risetimes{cc,ii}(:,1),falltimes{cc,ii}(:,1), rapidfreq{cc,ii}(:,1), trackLength_rapid{cc,ii}(:,1), trackLength_tot{cc,ii}(:,1), trackLength_ActiveFrac{cc,ii}(:,1), peaknums{cc,ii}(:,1), norm_peaknums{cc,ii}(:,1))')
end

end