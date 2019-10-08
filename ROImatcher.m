function varargout = ROImatcher(ROImatchCyt, ROImatchNuc, ...
     intDenWhole_ch1,intDenWhole_ch2,intDenWhole_ch3,intDenWhole_ch4, ...
       intDenNuc_ch1,  intDenNuc_ch2,  intDenNuc_ch3,  intDenNuc_ch4,...
       areaWhole_ch1,  areaWhole_ch2,  areaWhole_ch3,  areaWhole_ch4,...
         areaNuc_ch1,    areaNuc_ch2,    areaNuc_ch3,    areaNuc_ch4)

%Match CytROIs to list of NucROIs
[ROImatchNuc, ROImatchCyt] = cellpad(ROImatchNuc, ROImatchCyt);
[lgi idx] = cellfun(@(x,y) ismember(x,y), ROImatchNuc, ROImatchCyt, 'UniformOutput', 0);

%Rewrite all variables:
%1. Pad vectors to same length
%2. Delete NucROIs that don't exist in CytROIs
%3. Reorder CytROIs to match rows with NucROIs

[intDenNuc_ch1, intDenWhole_ch1] = cellpad(intDenNuc_ch1, intDenWhole_ch1);
intDenNuc_ch1   = cellfun(@(x,y)   x(y>0),    intDenNuc_ch1,   lgi, 'UniformOutput',0 );
intDenWhole_ch1 = cellfun(@(x,y,z) x(z(y>0)), intDenWhole_ch1, lgi, idx, 'UniformOutput',0);

[intDenNuc_ch2, intDenWhole_ch2] = cellpad(intDenNuc_ch2, intDenWhole_ch2);
intDenNuc_ch2   = cellfun(@(x,y)   x(y>0),    intDenNuc_ch2,   lgi, 'UniformOutput',0 );
intDenWhole_ch2 = cellfun(@(x,y,z) x(z(y>0)), intDenWhole_ch2, lgi, idx, 'UniformOutput',0);

[intDenNuc_ch3, intDenWhole_ch3] = cellpad(intDenNuc_ch3, intDenWhole_ch3);
intDenNuc_ch3   = cellfun(@(x,y)   x(y>0),    intDenNuc_ch3,   lgi, 'UniformOutput',0 );
intDenWhole_ch3 = cellfun(@(x,y,z) x(z(y>0)), intDenWhole_ch3, lgi, idx, 'UniformOutput',0);

[intDenNuc_ch4, intDenWhole_ch4] = cellpad(intDenNuc_ch4, intDenWhole_ch4);
intDenNuc_ch4   = cellfun(@(x,y)   x(y>0),    intDenNuc_ch4,   lgi, 'UniformOutput',0 );
intDenWhole_ch4 = cellfun(@(x,y,z) x(z(y>0)), intDenWhole_ch4, lgi, idx, 'UniformOutput',0);

[areaNuc_ch1, areaWhole_ch1] = cellpad(areaNuc_ch1, areaWhole_ch1);
areaNuc_ch1   = cellfun(@(x,y)   x(y>0),    areaNuc_ch1,   lgi, 'UniformOutput',0 );
areaWhole_ch1 = cellfun(@(x,y,z) x(z(y>0)), areaWhole_ch1, lgi, idx, 'UniformOutput',0);

[areaNuc_ch2, areaWhole_ch2] = cellpad(areaNuc_ch2, areaWhole_ch2);
areaNuc_ch2   = cellfun(@(x,y)   x(y>0),    areaNuc_ch2,   lgi, 'UniformOutput',0 );
areaWhole_ch2 = cellfun(@(x,y,z) x(z(y>0)), areaWhole_ch2, lgi, idx, 'UniformOutput',0);

[areaNuc_ch3, areaWhole_ch3] = cellpad(areaNuc_ch3, areaWhole_ch3);
areaNuc_ch3   = cellfun(@(x,y)   x(y>0),    areaNuc_ch3,   lgi, 'UniformOutput',0 );
areaWhole_ch3 = cellfun(@(x,y,z) x(z(y>0)), areaWhole_ch3, lgi, idx, 'UniformOutput',0);

[areaNuc_ch4, areaWhole_ch4] = cellpad(areaNuc_ch4, areaWhole_ch4);
areaNuc_ch4   = cellfun(@(x,y)   x(y>0),    areaNuc_ch4,   lgi, 'UniformOutput',0 );
areaWhole_ch4 = cellfun(@(x,y,z) x(z(y>0)), areaWhole_ch4, lgi, idx, 'UniformOutput',0);

%Outputs
varargout{1}  = intDenWhole_ch1;
varargout{2}  = intDenWhole_ch2;
varargout{3}  = intDenWhole_ch3;
varargout{4}  = intDenWhole_ch4;
varargout{5}  = intDenNuc_ch1;
varargout{6}  = intDenNuc_ch2;
varargout{7}  = intDenNuc_ch3;
varargout{8}  = intDenNuc_ch4;
varargout{9}  = areaWhole_ch1;
varargout{10} = areaWhole_ch2;
varargout{11} = areaWhole_ch3;
varargout{12} = areaWhole_ch4;
varargout{13} = areaNuc_ch1;
varargout{14} = areaNuc_ch2;
varargout{15} = areaNuc_ch3;
varargout{16} = areaNuc_ch4;

end