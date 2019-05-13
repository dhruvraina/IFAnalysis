%Plots stacked bars for fate decision experiments including error bars

%User Inputs
inputTxtFile         = '/Users/draina/Desktop/E11_Communifate/Figures/C6_F4_InductionTimeFateRatios_ConsolidatedMatlab.txt';
Consolidated         = readtable(inputTxtFile, 'Delimiter', 'tab','Format', '%s%f%f%f%f%f%s' );
file.outpath         = '/Users/draina/Desktop/E11_Communifate/Figures';
file.slashtype       = '/';
file.experimentName  = 'C6 F4 ind times - consolidated';
plotflag.imageFormat = 'svg';

%Main:
file.codeparent         = cd;
pfinder                 = strfind(file.codeparent, file.slashtype);
file.codeparent         = file.codeparent(1:pfinder(end));
addpath([file.codeparent file.slashtype 'ExchangeTools']);

ExperimentList = unique(Consolidated.Experiment, 'stable');
TreatmentList  = unique(Consolidated.Treatment,  'stable');


%Averaging and Gathering Errors:
for ii = 1:length(TreatmentList)
    
    idx                   = strncmp(TreatmentList{ii}, Consolidated.Treatment, 3);
    avg.DoublePos(ii)     = mean(Consolidated.DoublePositive(idx));
    avg.DoubleNeg(ii)     = mean(Consolidated.DoubleNegative(idx));
    avg.SingleG6P(ii)     = mean(Consolidated.OnlyGATA6Pos__ExcludesDP_(idx));
    avg.SingleNaP(ii)     = mean(Consolidated.OnlyNANOGPos__ExcludesDP_(idx));
    avg.TreatmentName(ii) = TreatmentList(ii);
    
    ste.DoublePos(ii)     = std(Consolidated.DoublePositive(idx))           /sqrt(length(ExperimentList));
    ste.DoubleNeg(ii)     = std(Consolidated.DoubleNegative(idx))           /sqrt(length(ExperimentList));
    ste.SingleG6P(ii)     = std(Consolidated.OnlyGATA6Pos__ExcludesDP_(idx))/sqrt(length(ExperimentList));
    ste.SingleNaP(ii)     = std(Consolidated.OnlyNANOGPos__ExcludesDP_(idx))/sqrt(length(ExperimentList));
    clear idx
end

resvec = vertcat(avg.SingleG6P,avg.DoublePos, avg.SingleNaP, avg.DoubleNeg);
stderr = vertcat(ste.SingleG6P,ste.DoublePos, ste.SingleNaP, ste.DoubleNeg);
errXid = repmat([1:length(TreatmentList)], length(ExperimentList),1);       %#ok<NBRAK>

%Change orientation:
resvec = resvec';
stderr = stderr';
errXid = errXid';

%Plotting:
fig2 = figure
h    = bar(resvec, 0.9, 'stacked');
hold on

%Set Properties
set(gca, 'XTick', 1:length(TreatmentList), 'XTickLabels', TreatmentList)
barTall = 10;
barWide = length(TreatmentList)*2;
xName   = '';
yName   = 'Fraction';

% Print text labels inside the bar graph
for i=1:size(resvec,1)
    for j=1:size(resvec,2)
        if resvec(i,j)>0.06                                          %Don't print values less than x% due to label space constraints
            labels_stacked = num2str((floor(resvec(i,j)*100))/100);  %Round digits to nearest lower full percent
            hText          = text(i, sum(resvec(i,1:j),2), labels_stacked);
            set(hText,'VerticalAlignment', 'top', 'HorizontalAlignment', 'left','FontSize',14, 'Color','w');
        end
    end
end

%Setting custom colours:
colmap = [136  40 144;  %Purple
           60 191 189;  %Blue
            0 126  61;  %Green
          241  90  41]; %Orange
colmap = colmap./256;   %convert from rgb space to [0 1]

h(1).FaceColor = colmap(1,:);
h(2).FaceColor = colmap(2,:);
h(3).FaceColor = colmap(3,:);
h(4).FaceColor = colmap(4,:);

%Errbar from FEX
errbar(errXid, cumsum(resvec,2),stderr, zeros(size(stderr)), 'w','linewidth', 5);

%Save + additional props:
xlabel(xName)
ylabel(yName)
title([file.experimentName '(std. Err)'])

fig2.PaperUnits    = 'inches';
fig2.PaperPosition = [0 0 barWide barTall];
set(gca,'FontSize', 14)
if ~isdir([file.outpath file.slashtype 'StackedBar'])
    mkdir([file.outpath file.slashtype 'StackedBar']);
end

switch plotflag.imageFormat
    case 'svg'
        print(fig2,[file.outpath file.slashtype 'StackedBar' file.slashtype 'Total'], '-painters', '-dsvg','-r200')
    case 'png'
        print(fig2,[file.outpath file.slashtype 'StackedBar' file.slashtype 'Total'], '-painters', '-dpng','-r200')
end




    
