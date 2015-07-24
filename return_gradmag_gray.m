function gradmag = return_gradmag_gray(BF)

hy = fspecial('sobel'); % filter for emphasizing vertical lines
hx = hy';               % filter for emphasizing horizontal lines
Iy = imfilter(double(BF), hy, 'replicate');
Ix = imfilter(double(BF), hx, 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);

