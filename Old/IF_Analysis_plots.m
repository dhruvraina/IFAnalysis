%2dsegalazyer
clearvars;
close all;

addpath('/Users/draina/Documents/MATLAB/Dhruv''s Stuff/MPI_Progs/ExchangeTools/raacampbell-notBoxPlot-2fbf98c/code');
addpath('/Users/draina/Documents/MATLAB/Dhruv''s Stuff/MPI_Progs/ExchangeTools/DrosteEffect-BrewerMap-a77e675')
addpath('/Users/draina/Documents/MATLAB/Dhruv''s Stuff/MPI_Progs/IFAnalysis');

slashtype = '/';


%% User Inputs:
mainDir = '/Users/draina/Desktop/kSpry4_NanogOct';
pathlist = {
        '/Users/draina/Desktop/kSpry4_NanogOct/AZ10'
        '/Users/draina/Desktop/kSpry4_NanogOct/AZ1000'
        '/Users/draina/Desktop/kSpry4_NanogOct/PD10'
        '/Users/draina/Desktop/kSpry4_NanogOct/PD100'
        '/Users/draina/Desktop/kSpry4_NanogOct/PD1000'
        '/Users/draina/Desktop/kSpry4_NanogOct/NT'
    
%     '/Users/draina/Desktop/2017_E2_2_eKTR_Kinetics/Single Slice Analysis/0Min'
%     '/Users/draina/Desktop/2017_E2_2_eKTR_Kinetics/Single Slice Analysis/1Min'
%     '/Users/draina/Desktop/2017_E2_2_eKTR_Kinetics/Single Slice Analysis/2Min'
%     '/Users/draina/Desktop/2017_E2_2_eKTR_Kinetics/Single Slice Analysis/3Min'
%     '/Users/draina/Desktop/2017_E2_2_eKTR_Kinetics/Single Slice Analysis/4Min'
%     '/Users/draina/Desktop/2017_E2_2_eKTR_Kinetics/Single Slice Analysis/8Min'
%     '/Users/draina/Desktop/2017_E2_2_eKTR_Kinetics/Single Slice Analysis/16Min'
  
% '/Users/draina/Desktop/2017_eKTR_ERK_Analysis/PDO3_Release'
% '/Users/draina/Desktop/2017_eKTR_ERK_Analysis/PDO3'
% '/Users/draina/Desktop/2017_eKTR_ERK_Analysis/ESL'

% '/Users/draina/Desktop/2017_June_ReleaseKinetics/1minRelease'
% '/Users/draina/Desktop/2017_June_ReleaseKinetics/2minRelease'
% '/Users/draina/Desktop/2017_June_ReleaseKinetics/3_5minRelease'
% '/Users/draina/Desktop/2017_June_ReleaseKinetics/7_5minRelease'
% '/Users/draina/Desktop/2017_June_ReleaseKinetics/15minRelease'
% '/Users/draina/Desktop/2017_June_ReleaseKinetics/30minRelease'
% % '/Users/draina/Desktop/2017_June_ReleaseKinetics/ESL'
% '/Users/draina/Desktop/2017_June_ReleaseKinetics/PDO3'

%'/Users/draina/Desktop/2017_eKTR_ERK/ESL'



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
% %    'ESL'
%      'PDO3'

% 'PDO3 Release (Max. ppERK)'
% 'PDO3 (Min. ppERK)'
% 'Serum+LIF'

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
        'NANOG'
        'VENUS'
        'OCT3/4'
    
%     'ppERK'
%     'Sensor'
%     'DAPI'
%     'Trans'
%     
%     'DAPI'
%     'Trans'
%     'Sensor'
%     'ppErk'
% %  


    
    };

% Reorder Channels according to how you want to plot them (expects 4
% channels, so just duplicate or set one to garbage channel)
ReorderChan = [ 3; 2; 1; 4 ];
%ReorderChan = [ 1; 3; 4; 2 ];

plotxlabel = [ChannelLabel(ReorderChan)];

save()

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
    
    %Copy totMed to the Normalization variable:
    totMedNorm_2 = totMed_2;
    totMedNorm_3 = totMed_3;
    totMedNorm_4 = totMed_4;
    
    % - 24.07.17: Why are you doing this? It's symmetrizing the distribution! Essentially this will tell you the same thing as S.D/Var only *less* quantitatively! REJECT.
    %Normalize to Mean:
    totMedNorm_2(:,1) = totMed_2(:,1)./mean(totMed_2(:,1));
    totMedNorm_3(:,1) = totMed_3(:,1)./mean(totMed_3(:,1));
    totMedNorm_4(:,1) = totMed_4(:,1)./mean(totMed_4(:,1));
    %%
    
    %Concatenate the results you want to see:
    resvec_norm{ctr2} = [totMedNorm_2; totMedNorm_3; totMedNorm_4];
    resvec_raw{ctr2}  = [totMed_2; totMed_3; totMed_4];
    
    %Per Treatment Box Plots (i.e. across all channels for one treatment group)
    ff = figure
    notBoxPlot(resvec_raw{ctr2}(:,1), resvec_raw{ctr2}(:,2))
    set(gca, 'xTickLabel', plotxlabel(2:4))
    ylim([0 100])
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
resscat = {};
colors = brewermap(length(resvec_raw), 'Dark2');

%Pull Out results for ALL TREATMENTS channel-by-channel
for ctr6 = 1:length(resvec_raw)
    
    %Extracting scatter points:
    resscat{ctr6,1} = resvec_raw{ctr6}(resvec_raw{ctr6}(:,2)==scatxCh,1);
    resscat{ctr6,2} = resvec_raw{ctr6}(resvec_raw{ctr6}(:,2)==scatyCh,1);
    
    %Plotting Scatter:
    qq = figure
    scatter(resscat{ctr6,1},resscat{ctr6,2}, ...
            50,'filled',...
            'MarkerFaceColor', colors(ctr6,:),...
            'MarkerFaceAlpha',3/9)
        
    % Graph Formatting and saving:
    xlim([0 90])
    ylim([0 90])
    xlabel([ChannelLabel{ReorderChan(scatxCh)} ' intensity (a.u.)'])
    ylabel([ChannelLabel{ReorderChan(scatyCh)} ' intensity (a.u.)'])
    title(pathlist_labels{ctr6})
    legend([pathlist_labels{ctr6} ': Cells Recorded: ' num2str(length(resscat{ctr6,1}))])
    
    if ~isdir([file(ctr6).maindir slashtype 'Scatter'])
        mkdir([file(ctr6).maindir slashtype 'Scatter']);
    end
    
    qq.PaperUnits = 'inches';
    qq.PaperPosition = [0 0 10 8];
    print(qq,[file(ctr6).maindir slashtype 'Scatter' slashtype 'sca_' file(ctr6).name ], '-painters', '-dsvg','-r200')
    close gcf
     
    
    
%     %Normalizing Scatter points:   
%     resscatNorm{ctr6,1} = resscat{ctr6,1}./median(resscat{ctr6,1});
%     resscatNorm{ctr6,2} = resscat{ctr6,2}./median(resscat{ctr6,2});
%     
%     %Plotting Normalized Scatter
%     gg = figure
%     scatter(resscatNorm{ctr6,1},resscatNorm{ctr6,2}, 'filled', 'MarkerFaceColor', colors(ctr6,:))
%     %xlim([0 2.5])
%     %ylim([0 2.5])
%     xlabel(['Median-Normalized ' ChannelLabel{ReorderChan(scatxCh)} ' intensity'])
%     ylabel(['Median-Normalized ' ChannelLabel{ReorderChan(scatyCh)} ' intensity'])
%     title(pathlist_labels{ctr6})
%     if ~isdir([file(ctr6).maindir slashtype 'NormalizedScatter'])
%         mkdir([file(ctr6).maindir slashtype 'NormalizedScatter']);
%     end
%     print(gg,[file(ctr6).maindir slashtype 'NormalizedScatter' slashtype 'Nsc_' file(ctr6).name], '-painters', '-dsvg','-r200')
%     close gcf
%     
     
end


%% Consolidated scatter:
scatxCh = 2;
scatyCh = 3;
resscat = {};
colors = brewermap(length(resvec_raw), 'Dark2');
oo = figure
for ctr6 = 1:length(resvec_raw)
    resscat{ctr6,1} = resvec_raw{ctr6}(resvec_raw{ctr6}(:,2)==scatxCh,1);   %x axis
    resscat{ctr6,2} = resvec_raw{ctr6}(resvec_raw{ctr6}(:,2)==scatyCh,1);   %y axis
    
    %Normalizing values - symmetrizes distributions. This is useless. DELETE!!:
    resscatNorm{ctr6,1} = resscat{ctr6,1}./median(resscat{ctr6,1});
    resscatNorm{ctr6,2} = resscat{ctr6,2}./median(resscat{ctr6,2});
    
    %plotting Scatter
    scatter(resscat{ctr6,1},resscat{ctr6,2}, ...
            50,'filled',...
            'MarkerFaceColor', colors(ctr6,:),...
            'MarkerFaceAlpha',3/9)
   xlim([0 80])
   ylim([0 90])   
    hold on
    legendary{ctr6} = [pathlist_labels{ctr6} ': Cells Recorded: ' num2str(length(resscat{ctr6,1}))];
   % keyboard
 end


legend(legendary)
if ~isdir([file(ctr6).maindir ])
    mkdir([file(ctr6).maindir ]);
end
xlabel(['Median Nuclear Intensity: ' ChannelLabel{ReorderChan(scatxCh)} ' channel (a.u.)'])
ylabel(['Median Nuclear Intensity: ' ChannelLabel{ReorderChan(scatyCh)} ' channel (a.u.)'])
    oo.PaperUnits = 'inches';
    oo.PaperPosition = [0 0 10 8];

print(oo,[file(1).maindir slashtype 'Scatter_Total'], '-painters', '-dsvg','-r200')
print(oo,[file(1).maindir slashtype 'Scatter_Total'], '-painters', '-dpng','-r200')


    