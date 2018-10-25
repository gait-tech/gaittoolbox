buf1 = quatmultiply(axang2quat([0 0 1 deg2rad(DEGRANGE(1))]), obj.(n).ori);
buf2 = quatmultiply(quatconj(obj.Pelvis.ori), buf1);
buf3 = quatrotate(quatconj(buf2), obj.(n).gyr);

theta_y = buf3(:,2);
err1(j) = skewness(theta_y);
%             err2(j) = peak2peak(theta_x);
updateFigureContents('Theta');
plot(theta_y)

theta_y0 = mean(theta_y(1:300));

[pks1, locs1] = findpeaks(theta_y, 'MinPeakHeight', mean(theta_y(1:300))+deg2rad(5));
[pks2, locs2] = findpeaks(-theta_y, 'MinPeakHeight', mean(-theta_y(1:300))+deg2rad(5));
% [pks1, locs1] = findpeaks(theta_y, 'MinPeakHeight', 0.1*theta_y0 + 0.1*max(theta_y));
% [pks2, locs2] = findpeaks(-theta_y, 'MinPeakHeight', -0.1*theta_y0 + 0.1*max(-theta_y));

updateFigureContents('Error');
plot(DEGRANGE, err1);


buf1 = quatmultiply(axang2quat([0 0 1 deg2rad(DEGRANGE(2))]), ori.(n));
buf2 = quatmultiply(quatconj(ori.Pelvis), buf1);
[theta_y theta_x theta_z] = quat2angle(buf2, 'YXZ');
%             err(j) = sum(theta_y.^2);
err1(j) = peak2peak(theta_y);
%             err2(j) = peak2peak(theta_x);

theta_y0 = mean(theta_y(1:500));
[pks1, locs1] = findpeaks(theta_y, 'MinPeakHeight', theta_y0+deg2rad(5));
[pks2, locs2] = findpeaks(-theta_y, 'MinPeakHeight', theta_y0+deg2rad(5));
            
updateFigureContents('Theta');
plot(rad2deg(theta_y))

updateFigureContents('Error');
plot(DEGRANGE, err1, DEGRANGE, err2);

up