%Plots stacked bars for fate decision experiments including error bars

clear all; close all;                                                      %#ok<CLALL>

%User Inputs
inputTxtFile         = '/Users/draina/Desktop/E11_Communifate/InductionTimes/Consolidated - logThresholds/C6_FateRatios_logThreshConsolidated.txt'; 
Consolidated         = readtable(inputTxtFile, 'Delimiter', 'tab','Format', '%s%f%f%f%f%f%s' );
file.outpath         = '/Users/draina/Desktop/E11_Communifate/InductionTimes/Consolidated - logThresholds/';
file.slashtype       = '/';
file.experimentName  = '191103logThresh';
plotflag.imageFormat = 'svg';
printNumbers         = 1; %Print values inside graph
minStringLenMatch    = 4; %Match first 4 characters from the Treatment Name
barTall              = 10; 
textSize             = 22;
confidenceInterval   = 1.96;  %1.96*sigma for 95%, 2*sigma for 95.45%
ylimAuto             = 1;
ylim1                = 0;
ylim2                = 100;

%Main:
file.codeparent         = cd;
pfinder                 = strfind(file.codeparent, file.slashtype);
file.codeparent         = file.codeparent(1:pfinder(end));
addpath([file.codeparent file.slashtype 'ExchangeTools']);

ExperimentList = unique(Consolidated.Experiment, 'stable');
TreatmentList  = unique(Consolidated.Treatment,  'stable');

%Find shortest name in TreatmentList for strmatch operations:
%[minStringLen, id] = min(cellfun(@(x) length(x), TreatmentList));

%Averaging and Gathering Errors:
for ii = 1:length(TreatmentList)
    
    idx                   = strncmp(TreatmentList{ii}, Consolidated.Treatment, minStringLenMatch);
    avg.DoublePos(ii)     = mean(Consolidated.DoublePositive(idx));
    avg.DoubleNeg(ii)     = mean(Consolidated.DoubleNegative(idx));
    avg.SingleG6P(ii)     = mean(Consolidated.OnlyGATA6Pos__ExcludesDP_(idx));
    avg.SingleNaP(ii)     = mean(Consolidated.OnlyNANOGPos__ExcludesDP_(idx));
    avg.TreatmentName(ii) = TreatmentList(ii);
    
    ste.DoublePos(ii)     = confidenceInterval * std(Consolidated.DoublePositive(idx))           /sqrt(length(ExperimentList));
    ste.DoubleNeg(ii)     = confidenceInterval * std(Consolidated.DoubleNegative(idx))           /sqrt(length(ExperimentList));
    ste.SingleG6P(ii)     = confidenceInterval * std(Consolidated.OnlyGATA6Pos__ExcludesDP_(idx))/sqrt(length(ExperimentList));
    ste.SingleNaP(ii)     = confidenceInterval * std(Consolidated.OnlyNANOGPos__ExcludesDP_(idx))/sqrt(length(ExperimentList));
    clear idx
end

resvec = vertcat(avg.SingleG6P,avg.DoublePos, avg.SingleNaP, avg.DoubleNeg);
stderr = vertcat(ste.SingleG6P,ste.DoublePos, ste.SingleNaP, ste.DoubleNeg);
errXid = repmat([1:length(TreatmentList)], 4,1);       %#ok<NBRAK>

%Change orientation:
resvec = resvec';
stderr = stderr';
errXid = errXid';

%Plotting:
fig2 = figure
h    = bar(resvec, 0.75, 'stacked');
hold on

%Set Properties
set(gca, 'XTick', 1:length(TreatmentList), 'XTickLabels', TreatmentList)
%barTall = 20;
barWide = length(TreatmentList)*2;
xName   = '';
yName   = 'Fraction';

% Print text labels inside the bar graph
if printNumbers
for i=1:size(resvec,1)                                                      %#ok<UNRCH>
    for j=1:size(resvec,2)
        if resvec(i,j)>0.06                                          %Don't print values less than x% due to label space constraints
            labels_stacked = num2str((floor(resvec(i,j)*100))/100);  %Round digits to nearest lower full percent
            hText          = text(i, sum(resvec(i,1:j),2), labels_stacked);
            set(hText,'VerticalAlignment', 'top', 'HorizontalAlignment', 'left','FontSize',14, 'Color','w');
        end
    end
end
end

%Setting custom colours:
colmap = [217 111 171;  %Pink     GATA
          210 188  44;  %Yellow   DoublePositive
            0 168 110;  %Green    NANOG
          140 198 236]; %Blue     DoubleNegative
colmap = colmap./256;   %convert from rgb space to [0 1]

h(1).FaceColor = colmap(1,:);
h(2).FaceColor = colmap(2,:);
h(3).FaceColor = colmap(3,:);
h(4).FaceColor = colmap(4,:);

%Errbar from FEX - for < R2015b
%errbar(errXid, cumsum(resvec,2),stderr, zeros(size(stderr)), 'b+-','linewidth', 1);

errorbar(errXid, cumsum(resvec,2), stderr, zeros(size(stderr)),'.k');
%Save + additional props:
xlabel(xName)
ylabel(yName)
title([file.experimentName '(std. Err)'])

if ~ylimAuto;       ylim([ylim1 ylim2]);           end

fig2.PaperUnits    = 'inches';
fig2.PaperPosition = [0 0 barWide barTall];
set(gca,'FontSize', textSize)
if ~isdir([file.outpath file.slashtype 'StackedBar'])
    mkdir([file.outpath file.slashtype 'StackedBar']);
end

switch plotflag.imageFormat
    case 'svg'
        print(fig2,[file.outpath file.slashtype 'StackedBar' file.slashtype file.experimentName], '-painters', '-dsvg','-r200')
    case 'png'
        print(fig2,[file.outpath file.slashtype 'StackedBar' file.slashtype file.experimentName], '-painters', '-dpng','-r200')
end




    
