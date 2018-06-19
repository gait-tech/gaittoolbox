idx = 1:data.dataV.nSamples;
actBody = data.dataV.togrBody(idx, {});

updateFigureContents('Animation');
actBodyLimits = [actBody.xlim() actBody.ylim() actBody.zlim()]*3;
i = 1;
while i <= actBody.nSamples
    clf; grid;
    xlim(actBodyLimits(1:2)); 
    ylim(actBodyLimits(3:4)); 
    zlim(actBodyLimits(5:6));  
    xlabel('x'); ylabel('y'); zlabel('z');
    view(20, 30);
    pelib.viz.plotLowerBody(actBody, i, true, false);
    i = i+5;
    pause(1/1000);
end