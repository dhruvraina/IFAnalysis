%Plots stacked bars for fate decision experiments including error bars
clear all; close all;                                                      %#ok<CLALL>
addpath('/Users/draina/Documents/Code/MATLAB/DhruvTools/');


%User Inputs
inputTxtFile         = '/Users/draina/Desktop/E11_Communifate/DifferentClones/Consolidated/logThresh/AllExperiments_InclFGF4.txt';
Consolidated         = readtable(inputTxtFile, 'Delimiter', 'tab','Format', '%s%f%f%f%f%f%s' );
file.outpath         = '/Users/draina/Desktop/E11_Communifate/DifferentClones/Consolidated/logThresh';
file.slashtype       = '/';
file.experimentName  = '191212_ReDo';
plotflag.imageFormat = 'svg';
printNumbers         = 0;     %Print values inside graph
minStringLenMatch    = 5;     %Match first x characters from the Treatment Name
barTall              = 10;    %Height of plot
confidenceInterval   = 1.96;  %1.96*sigma for 95%, 2*sigma  for 95.45%
ylimAuto             = 1;     %Overrides ylims below if ON
ylim1                = 0;
ylim2                = 100;
textSize             = 12;

%Main:
file.codeparent         = cd;
pfinder                 = strfind(file.codeparent, file.slashtype);
file.codeparent         = file.codeparent(1:pfinder(end));
addpath([file.codeparent file.slashtype 'ExchangeTools']);

ExperimentList = unique(Consolidated.Experiment, 'stable');
TreatmentList  = unique(Consolidated.Treatment,  'stable');
nTreatments    = length(TreatmentList);


%Find shortest name in TreatmentList for strmatch operations:
%[minStringLen, id] = min(cellfun(@(x) length(x), TreatmentList));

%Averaging and Gathering Errors:
for ii = 1:nTreatments
    
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
barWide = length(TreatmentList)*2; %Scale plot width
xName   = '';
yName   = 'Fraction';

% Print text labels inside the bar graph
if printNumbers
    for i=1:size(resvec,1)                                                 %#ok<UNRCH>
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
errbar(errXid, cumsum(resvec,2),stderr, zeros(size(stderr)), 'k+-','linewidth', 1);
%errorbar(errXid, cumsum(resvec,2), stderr, zeros(size(stderr)),'.k');

%Save + additional props:
xlabel(xName)
ylabel(yName)
title([cleanNames({file.experimentName},'_') '(std. Err)'])

if ~ylimAuto;       ylim([ylim1 ylim2]);           end

fig2.PaperUnits    = 'inches';
fig2.PaperPosition = [0 0 barWide barTall];
set(gca,'FontSize', textSize)
set(gca,'TickDir', 'out')
if ~isdir([file.outpath file.slashtype 'StackedBar'])
    mkdir([file.outpath file.slashtype 'StackedBar']);
end


%Write means to file:
fileID = fopen([file.outpath file.slashtype file.experimentName '_Reformat.txt'],'w');
fprintf(fileID, '%1$s\t %2$s\t %3$s\t %4$s\t %5$s\r\n', 'Treatment', 'Mean GATA6Pos', 'Mean DP', 'Mean NANOGPos', 'Mean DN');

for tnum = 1:nTreatments
    fprintf(fileID,'%1$s\t', TreatmentList{tnum});
    fprintf(fileID,'%4.3f\t %4.3f\t %4.3f\t %4.3f\r\n', resvec(tnum,1), resvec(tnum,2), resvec(tnum,3), resvec(tnum,4));
end

%Write st.err to file:
fileID = fopen([file.outpath file.slashtype file.experimentName '_stErr.txt'],'w');
fprintf(fileID, '%1$s\t %2$s\t %3$s\t %4$s\t %5$s\r\n', 'Treatment', 'stdErr GATA6Pos', 'stdErr DP', 'stdErr NANOGPos', 'stdErr DN');

for tnum = 1:nTreatments
    fprintf(fileID,'%1$s\t', TreatmentList{tnum});
    fprintf(fileID,'%4.3f\t %4.3f\t %4.3f\t %4.3f\r\n', stderr(tnum,1), stderr(tnum,2), stderr(tnum,3), stderr(tnum,4));
end


%Print graphs
switch plotflag.imageFormat
    case 'svg'
        print(fig2,[file.outpath file.slashtype 'StackedBar' file.slashtype file.experimentName], '-painters', '-dsvg','-r200')
    case 'png'
        print(fig2,[file.outpath file.slashtype 'StackedBar' file.slashtype file.experimentName], '-painters', '-dpng','-r200')
end





