
function varargout = corrline(xv, yv)
[p,err] = polyfit(xv,yv,1);   % First order polynomial
y_fit = polyval(p,xv,err);   % Values on a line
y_dif = yv - y_fit;          % y value difference (residuals)
SSdif = sum(y_dif.^2);      % Sum square of difference
SStot = (length(yv)-1)*var(xv);   % Sum square of y taken from variance
rsq = 1-SSdif/SStot;
varargout{1} = xv;
varargout{2} = y_fit;
end

