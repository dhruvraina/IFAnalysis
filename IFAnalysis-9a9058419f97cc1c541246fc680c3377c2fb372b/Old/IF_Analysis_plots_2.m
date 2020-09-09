%2dsegalazyer
clearvars;
close all;

addpath('/Users/draina/Documents/MATLAB/dr_progs/MPI_Progs/ExchangeTools/raacampbell-notBoxPlot-2fbf98c/code');
addpath('/Users/draina/Documents/MATLAB/dr_progs/MPI_Progs/ExchangeTools/DrosteEffect-BrewerMap-a77e675')
addpath('/Users/draina/Documents/MATLAB/dr_progs/MPI_Progs/IFAnalysis');

slashtype = '/';

%% User Inputs:
mainDir = '/Users/draina/Desktop/Spry4_reporter_paper (CopiedToServer)/Raw Data (Copied to Server)/Immunofluorescence/2017_E6_kSpry4_NanogOct/Analysed Data';
loadmode = 1; %is on
if exist([mainDir slashtype 'IF_Analysis_Settings.mat']) && loadmode==1
    load([mainDir slashtype 'IF_Analysis_Settings.mat']);
else
    
    pathlist = {
        '/Users/draina/Desktop/Spry4_reporter_paper (CopiedToServer)/Raw Data (Copied to Server)/Immunofluorescence/2017_E6_kSpry4_NanogOct/Analysed Data/AZ10'
        '/Users/draina/Desktop/Spry4_reporter_paper (CopiedToServer)/Raw Data (Copied to Server)/Immunofluorescence/2017_E6_kSpry4_NanogOct/Analysed Data/AZ1000'
        '/Users/draina/Desktop/Spry4_reporter_paper (CopiedToServer)/Raw Data (Copied to Server)/Immunofluorescence/2017_E6_kSpry4_NanogOct/Analysed Data/PD10'
        '/Users/draina/Desktop/Spry4_reporter_paper (CopiedToServer)/Raw Data (Copied to Server)/Immunofluorescence/2017_E6_kSpry4_NanogOct/Analysed Data/PD100'
        '/Users/draina/Desktop/Spry4_reporter_paper (CopiedToServer)/Raw Data (Copied to Server)/Immunofluorescence/2017_E6_kSpry4_NanogOct/Analysed Data/PD1000'
        '/Users/draina/Desktop/Spry4_reporter_paper (CopiedToServer)/Raw Data (Copied to Server)/Immunofluorescence/2017_E6_kSpry4_NanogOct/Analysed Data/NT'
        };
    
    pathlist_labels = {
        'AZ10'
        'AZ1000'
        'PD10'
        'PD100'
        'PD1000'
        'NT'
        };
    
    
    %Order Channel Labels according to image
    ChannelLabel = {
        'DAPI'
        'Trans'
        'VENUS'
        'OCT4'
        'NANOG'
        };
    
    % Reorder Channels according to how you want to plot them (expects 4
    % channels, so just duplicate or set one to garbage channel)
    ReorderChan = [1 3 5 4];
    plotxlabel = [ChannelLabel(ReorderChan)];
    save([mainDir slashtype 'IF_Analysis_Settings.mat'], 'mainDir', 'pathlist', 'pathlist_labels', 'ReorderChan', 'plotxlabel', 'ChannelLabel', '-v7');
end
%% Main


%Reassign ChannelVec to easy to read vars
ch1 = ReorderChan(1);
ch2 = ReorderChan(2);
ch3 = ReorderChan(3);
ch4 = ReorderChan(4);


for ctr2 = 1:length(pathlist)
    file(ctr2).path    = pathlist{ctr2};
    file(ctr2).name    = file(ctr2).path(find(file(ctr2).path==slashtype, 1, 'last')+1:end);
    file(ctr2).maindir = file(ctr2).path(1:(find(file(ctr2).path==slashtype, 1, 'last')-1));
    
    %Initialize
    median_ch2 = {};
    median_ch3 = {};
    median_ch4 = {};
    
    %Only look for .txt files
    dir_cell = dir( fullfile(file(ctr2).path, '*.txt'));
    
    
    %readtable is faster than loaddata! Loading only Medians in here
    for c1 = 1:length(dir_cell)
        temptable      = readtable([file(ctr2).path slashtype dir_cell(c1).name], 'delimiter', '\t');
        median_ch2{c1} = temptable.Median(temptable.Ch==ch2);
        median_ch3{c1} = temptable.Median(temptable.Ch==ch3);
        median_ch4{c1} = temptable.Median(temptable.Ch==ch4);
    end
    
    
    %Concat. everything from loop
    totMed_2 = cell2mat(median_ch2');
    totMed_3 = cell2mat(median_ch3');
    totMed_4 = cell2mat(median_ch4');
    
    %Assign the Channel Indices
    totMed_2(1:length(totMed_2),2)=2;
    totMed_3(1:length(totMed_3),2)=3;
    totMed_4(1:length(totMed_4),2)=4;
    
    
    % Stuff here is for INDIVIDUAL 'all-channels-in-one-treatment' plots
    
    %Assign Plotting Indicies
    totMed_2(1:length(totMed_2),2)=2;  %2,3,4 are indicies that let me find ch2,3,4 later
    totMed_3(1:length(totMed_3),2)=3;
    totMed_4(1:length(totMed_4),2)=4;
    
    resvec_raw{ctr2}  = [totMed_2; totMed_3; totMed_4];
    
    %Per Treatment Box Plots (i.e. across all channels for one treatment group)
    ff = figure
    notBoxPlot(resvec_raw{ctr2}(:,1), resvec_raw{ctr2}(:,2))
    set(gca, 'xTickLabel', plotxlabel(2:4))
    ylim([0 250])
    ylabel('Median Intensities')
    title(file(ctr2).name)
    if ~isdir([file(ctr2).path slashtype 'Figs'])
        mkdir([file(ctr2).path slashtype 'Figs']);
    end
    print(ff,[file(ctr2).path slashtype 'Figs' slashtype 'Nbx_' file(ctr2).name ], '-painters', '-dsvg','-r200')
    close gcf
    
    
    %Clear Vars for loop:
    clear totMedNorm_2 totMedNorm_4 totMedNorm_3 totMed_2 totMed_3 totMed_4  temptable
end




%% Per-Channel Box Plots (i.e. across all treatment groups for one channel)
for ctr5 = 2:4
    analyzeCh = ReorderChan(ctr5);
    resfinal  = [];
    
    %Pull Out results -> For one channel at a time, show all treatments
    for ctr4 = 1:length(resvec_raw)
        tempvec                      = resvec_raw{ctr4}(resvec_raw{ctr4}(:,2)==ctr5,1);
        tempvec(1:length(tempvec),2) = ctr4;
        resfinal                     = [resfinal; tempvec];
        clear tempvec
    end
    resfinal(end+1,2) = 2;
    
    %Plotting
    tt = figure
    notBoxPlot(resfinal(:,1), resfinal(:,2))
    set(gca, 'XTickLabel', [pathlist_labels])
    title(['Median ' ChannelLabel{analyzeCh} ' Per Cell'])
    ylabel('Median Intensity')
    if ~isdir([file(1).maindir slashtype 'BoxPlot'])
        mkdir([file(1).maindir slashtype 'BoxPlot']);
    end
    print(tt,[file(1).maindir slashtype 'BoxPlot' slashtype 'box_' ChannelLabel{analyzeCh}], '-painters', '-dsvg','-r200')
    close gcf
    
    
end





%% Plot Scatter, which channel vs. which channel?
scatxCh = 2;
scatyCh = 3;
invertX = 0; %inverts the X axis values
invertY = 0; %inverts the Y axis values
corrprint = 1;
resscat = {};
colors = brewermap(length(resvec_raw), 'Dark2');
q =1;
%Pull Out results for ALL TREATMENTS channel-by-channel
for ctr6 = [6 3 4 5] %1:length(resvec_raw) % %
    
    %Extracting scatter points:
    if invertX ==1
        resscat{ctr6,1} = 1./(resvec_raw{ctr6}(resvec_raw{ctr6}(:,2)==scatxCh,1));
    else
        resscat{ctr6,1} = resvec_raw{ctr6}(resvec_raw{ctr6}(:,2)==scatxCh,1);
    end
    
    if invertY ==1
        resscat{ctr6,2} = 1./(resvec_raw{ctr6}(resvec_raw{ctr6}(:,2)==scatyCh,1));
    else
        resscat{ctr6,2} = resvec_raw{ctr6}(resvec_raw{ctr6}(:,2)==scatyCh,1);
    end
    
    resscat_mean(q,1)  =  mean(resscat{ctr6,1});
    resscat_mean(q,2)  =  mean(resscat{ctr6,2});
    resscat_xyerr(q,1) =  std(resscat{ctr6,1});
    resscat_xyerr(q,2) =  std(resscat{ctr6,2});
    %Plotting Scatter:
    qq = figure
    scatter(resscat{ctr6,1},resscat{ctr6,2}, ...
        50,'filled',...
        'MarkerFaceColor', colors(q,:),...
        'MarkerFaceAlpha',3/9)
    q = q+1;
    % Graph Formatting and saving:
     ylim([0 250])
     xlim([0 250])
    
    if invertX==1
        invXstring = '1/';
    else
        invXstring = '';
    end
    
    if invertY==1
        invYstring = '1/';
    else
        invYstring = '';
    end
    
    xlabel([invXstring ChannelLabel{ReorderChan(scatxCh)} ' intensity (a.u.)'])
    ylabel([invYstring ChannelLabel{ReorderChan(scatyCh)} ' intensity (a.u.)'])
    title(pathlist_labels{ctr6})
    legend([pathlist_labels{ctr6} ': Cells Recorded: ' num2str(length(resscat{ctr6,1}))])
    
    if corrprint ==1
        hold on
        [clinx cliny] = corrline(resscat{ctr6,1}, resscat{ctr6,2});
        [pearsn pval] = corr(resscat{ctr6,1}, resscat{ctr6,2});
        rsq = pearsn^2;
        plot(clinx, cliny)
        plotstrg = {['r^2 = ' num2str(rsq)], ['p-value = ' num2str(pval)]};
        dim = [.2 .5 .3 .3];
        annotation('textbox',dim,'String',plotstrg,'FitBoxToText','on');
    end
    
    if ~isdir([file(ctr6).maindir slashtype 'Scatter'])
        mkdir([file(ctr6).maindir slashtype 'Scatter']);
    end
    
    qq.PaperUnits = 'inches';
    qq.PaperPosition = [0 0 10 8];
    print(qq,[file(ctr6).maindir slashtype 'Scatter' slashtype 'sca_' file(ctr6).name ], '-painters', '-dsvg','-r200')
    close gcf
    
    
    clear clinx cliny
end


%% Consolidated scatter:
scatxCh = 2;
scatyCh = 3;
resscat = {};
colors = brewermap(length(resvec_raw), 'Dark2');
oo = figure
reorderTreatment = [6 3 4 5];
corrmat = [];
invertX = 0;
invertY = 0;
medianOnly = 0;

for ctr6 = 1:length(reorderTreatment)
    ww = reorderTreatment(ctr6);
    
    if invertX ==1 && medianOnly ==0
        resscat{ctr6,1} = 1./(resvec_raw{ww}(resvec_raw{ww}(:,2)==scatxCh,1));
    else
        resscat{ctr6,1} = resvec_raw{ww}(resvec_raw{ww}(:,2)==scatxCh,1);
    end
    
    if invertY ==1 && medianOnly ==0
        resscat{ctr6,2} = 1./(resvec_raw{ww}(resvec_raw{ww}(:,2)==scatyCh,1));   %y axis
    else
        resscat{ctr6,2} = resvec_raw{ww}(resvec_raw{ww}(:,2)==scatyCh,1);   %y axis
    end
    
    if medianOnly==1
        resscat{ctr6,1} = median(resvec_raw{ww}(resvec_raw{ww}(:,2)==scatxCh,1)); %x axis
        resscat{ctr6,2} = median((resvec_raw{ww}(resvec_raw{ww}(:,2)==scatyCh,1)));   %y axis
    end
    
    tempvec = resscat{ctr6,1}';
    tempvec = [tempvec; resscat{ctr6,2}'];
    corrmat = [corrmat; tempvec'];
    %plotting Scatter
    scatter(resscat{ctr6,1},resscat{ctr6,2}, ...
        50,'filled',...
        'MarkerFaceColor', colors(ctr6,:),...
        'MarkerFaceAlpha',3/9)
    hold on
    legendary{ctr6} = [pathlist_labels{ww} ': Cells Recorded: ' num2str(length(resscat{ctr6,1}))];
    
    clear tempvec
end
ylim([0 250])
xlim([0 250])


if corrprint ==1
    [clinx cliny] = corrline(corrmat(:,1), corrmat(:,2));
    [pearsn pval] = corr(corrmat(:,1), corrmat(:,2));
    rsq = pearsn^2;
    plot(clinx, cliny)
    plotstrg = {['r^2 = ' num2str(rsq)], ['p-value = ' num2str(pval)]};
    dim = [.2 .5 .3 .3];
    annotation('textbox',dim,'String',plotstrg,'FitBoxToText','on');
end


legend(legendary)
if ~isdir([file(ww).maindir ])
    mkdir([file(ww).maindir ]);
end

if invertX==1
    invXstring = '1/';
else
    invXstring = '';
end

if invertY==1
    invYstring = '1/';
else
    invYstring = '';
end

if medianOnly==1
    medstring = '_median';
else
    medstring = '';
end

xlabel([invXstring 'Median Nuclear Intensity: ' ChannelLabel{ReorderChan(scatxCh)} ' channel (a.u.)'])
ylabel([invYstring 'Median Nuclear Intensity: ' ChannelLabel{ReorderChan(scatyCh)} ' channel (a.u.)'])
oo.PaperUnits = 'inches';
oo.PaperPosition = [0 0 10 8];

print(oo,[file(1).maindir slashtype 'Scatter_Total_3' medstring], '-painters', '-dsvg','-r200')
print(oo,[file(1).maindir slashtype 'Scatter_Total_3' medstring], '-painters', '-dpng','-r200')

