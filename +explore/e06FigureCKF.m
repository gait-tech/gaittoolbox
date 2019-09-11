% ======================================================================
%> @brief Generate figure for algorithm overview
%> @author Luke Sy
% ======================================================================
body = pelib.grBody();

body.xyzColor = {'r', 'g', 'b'};
body.rplColor = {'k', 'k', 'k'};
body.lnSymbol = '-';
body.ptSymbol = 'o';
        
% body position
body.MIDPEL = [0 0 0.4];
body.LFEP = [0 0.2 0.4];
body.LFEO = [0.05 0.2 0.2];
body.LTIO = [0 0.2 0];
body.RFEP = [0 -0.2 0.4];
body.RFEO = [0.05 -0.2 0.2];
body.RTIO = [0 -0.2 0];

% body orientation
body.qRPV = rotm2quat(eye(3));
qLTH_y = [0 1 0];
qLTH_z = (body.LFEP - body.LFEO) / norm(body.LFEP - body.LFEO);
body.qLTH = rotm2quat([cross(qLTH_y, qLTH_z)' qLTH_y' qLTH_z']);
qLSK_z = (body.LFEO - body.LTIO) / norm(body.LFEO - body.LTIO);
body.qLSK = rotm2quat([cross(qLTH_y, qLSK_z)' qLTH_y' qLSK_z']);
qRTH_z = (body.RFEP - body.RFEO) / norm(body.RFEP - body.RFEO);
body.qRTH = rotm2quat([cross(qLTH_y, qRTH_z)' qLTH_y' qRTH_z']);
qRSK_z = (body.RFEO - body.RTIO) / norm(body.RFEO - body.RTIO);
body.qRSK = rotm2quat([cross(qLTH_y, qRSK_z)' qLTH_y' qRSK_z']);

% plot
clf; hold;
nOrig = repmat([0.5 -0.2 0], [3, 1]); nAxis = eye(3); nColor = 'rgb';
for i=1:3
    quiver3(nOrig(i,1), nOrig(i,2), nOrig(i,3), ...
            nAxis(i,1), nAxis(i,2), nAxis(i,3), ...
            0.15, 'Color', nColor(i), 'LineWidth', 2, 'MaxHeadSize', 1);
end
pelib.viz.plotLowerBody(body, 1, 2, false); view(132, 17);
xlim([-0.1 0.7]); ylim([-0.4 0.4]); zlim([-0.1 0.7]); axis square;

print(gcf, 'body-skeleton-orig', '-dtiff', '-r300');