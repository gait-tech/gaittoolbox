% Acc Exploration
fs=100;
estAccMP1 = estAcc(:,1:3); estAccLA1 = estAcc(:,4:6); estAccRA1 = estAcc(:,7:9); 
actAccMP1 = actAcc(:,1:3); actAccLA1 = actAcc(:,4:6); actAccRA1 = actAcc(:,7:9);

updateFigureContents('Hip Acc');
clf; grlib.viz.plotComparison(estAccMP1, actAccMP1, fs);

updateFigureContents('LA Acc');
clf; grlib.viz.plotComparison(estAccLA1, actAccLA1, fs);

updateFigureContents('RA Acc');
clf; grlib.viz.plotComparison(estAccRA1, actAccRA1, fs);
