load('C:\Users\z5151460\workspace\gaitrecon\eXsens\xsens-s1-mt_012000eb-004-Drx+M03+C222.mat')
updateFigureContents('Position1');
actBodyRel = actBody.changeRefFrame('MIDPEL');
actBodyRel.xyzColor = {'k', 'k', 'k'};
actBodyRel.lnSymbol = '--';
estBodyRel = estBody.changeRefFrame('MIDPEL');
estBodyRel.lnSymbol = '-';
grlib.viz.plotPosition({estBodyRel, actBodyRel}, {'RTIO'});

updateFigureContents('Animation');
xlabel('x'); ylabel('y'); zlabel('z');
estBodyRel = estBody.changeRefFrame('MIDPEL');
actBodyRel = actBody.changeRefFrame('MIDPEL');
actBodyRel.xyzColor = {'k', 'k', 'k'};
estBody2 = estBodyRel.toWorldFrame(actBody.MIDPEL, actBody.qRPV);
actBody2 = actBodyRel.toWorldFrame(actBody.MIDPEL+[1 0 0], actBody.qRPV);
estBodyLimits = [estBody2.xlim()+[-1 +1] estBody2.ylim()+[-1 +1] estBody2.zlim()];
i = 1; az = 0; el = 180;
while i <= estBody.nSamples
    [az, el] = view;
    clf; grid;
    xlim(estBodyLimits(1:2)); 
    ylim(estBodyLimits(3:4)); 
    zlim(estBodyLimits(5:6));  
    xlabel('x'); ylabel('y'); zlabel('z');
    view(az, el);
    grlib.viz.plotLowerBody(estBody2, i, true, false);
    grlib.viz.plotLowerBody(actBody2, i, true, false);
    i = i+10;
    pause(1/1000);
end

% crouch drift
load('C:\Users\z5151460\workspace\gaitrecon\eXsens\xsens-s1-mt_012000eb-004-Dsx+M03+C002.mat')

updateFigureContents('Animation');
xlabel('x'); ylabel('y'); zlabel('z');
estBodyRel = estBody.changeRefFrame('MIDPEL');
actBodyRel = actBody.changeRefFrame('MIDPEL');
estBody2 = estBodyRel.toWorldFrame(actBody.MIDPEL, actBody.qRPV);
actBody2 = actBodyRel.toWorldFrame(actBody.MIDPEL+[1 0 0], actBody.qRPV);
estBodyLimits = [estBody2.xlim()+[-1 +1] estBody2.ylim()+[-1 +1] estBody2.zlim()];
i = 1500; az = 0; el = 180;

[az, el] = view;
clf; grid;
xlim([3 5]); 
ylim([0.5 2.5]); 
zlim(estBodyLimits(5:6));  
xlabel('x'); ylabel('y'); zlabel('z');
view(az, el);
grlib.viz.plotLowerBody(estBody2, i, true, false);
grlib.viz.plotLowerBody(actBody2, i, true, false);
i = i+20;

    