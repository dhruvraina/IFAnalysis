%Interactive Nuc

function varargout = internuc(nuc, bw4)
nucsize = 3;

hy = fspecial('sobel');
hx = hy';
Iy = imfilter(double(nuc), hy, 'replicate');
Ix = imfilter(double(nuc), hx, 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);

se = strel('disk', nucsize);
Ie = imerode(gradmag, se);
Iobr = imreconstruct(Ie, nuc);
Iobrd = imdilate(Iobr, se);

nuc2 = Iobrd.*bw4;
gg = im2bw(nuc2, graythresh(nuc2));
che = imerode(gg, strel('disk', 1));
nucprop = regionprops(che, 'Area');
nucfilt = ceil(mean([nucprop(:).Area])/3);

imclear = bwareaopen(che, nucfilt);
[lbl num] = bwlabel(imclear);

%regmax
varargout{1} = lbl;
varargout{2} = num;