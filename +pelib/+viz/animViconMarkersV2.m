% ======================================================================
%> @brief quickly visualise the vicon markers
%> @author Michael Del Rosario (UNSW GSBME)
%>
%> Source: https://github.com/mikedelr/Vicon-2-MATLAB-tools/blob/master/animViconMarkersV2.m
%>
%> Example usage:
%>      animViconMarkers('markerData',markerData,'markerSet',markerSet)
%>
%> @param bodys Body instance(s) to be plotted
%> @param parts String(s) of body point(s) to be plotted.
%>
%> @param markerData struct with fields:
%>           .Names - a 'containers.Map' java object whose keys are the
%>                    marker names, and whose values denote the columns
%>                    in the (.Pos) matrix denote the co-ordinates in R3
%>                    of the markers (units in mm)
%>           .Pos   - a matrix containing the 3D co-ordinates of each
%>                    marker
%> @param markerSet cell array of size [n_rows, 2 columns]
%> @param markerSetColour cell array of colour markers, e.g., {'b';'r';'g'}
% ======================================================================
function [hFig ] = animViconMarkersV2( varargin )
bMarkerSet = false;
bMarkerColours = false;
markerSet  = [];
speedIncrement = 1;
bSpin = false;
globalIdx = []; 
figstr = 'Vicon Markers';
for i=1:2:nargin
    if  strcmp(varargin{i}, 'markerData'),
        markerData = varargin{i+1};
    elseif  strcmp(varargin{i}, 'markerSet'),
        markerSet = varargin{i+1};
        bMarkerSet = true;
    elseif strcmp(varargin{i}, 'markerSetColour')
        markerSetColour = varargin{i+1};
        bMarkerColours = true;
    elseif  strcmp(varargin{i}, 'figname'),
        figstr = strcat(figstr,'_',varargin{i+1});
    elseif strcmp(varargin{i}, 'speed'),
        speedIncrement = varargin{i+1};
    elseif strcmp(varargin{i}, 'spin'),
        bSpin = varargin{i+1};
    else error('Invalid argument');
    end    
end

txtOffset = 10;
hFig = figure('name',figstr,'Tag','MasterViconAnim');whitebg(hFig,[0.1 0.1 0.1])
az=-66;el=26;view(az,el);
set(hFig,'toolbar','figure')
DEFAULT_PLAYUNTILFRAMESKIP = 1000;

%% ------------ Determine suitable axis limits for 3D plot
[NUM_FRAMES,NCOLS]=size(markerData.Pos);
xyzIDX = 3:NCOLS;
NUM_MARKS = double(markerData.Names.Count);
XYZ = zeros(NUM_MARKS*NUM_FRAMES,3);
for m = 1:NUM_MARKS
    idx = bsxfun(@plus,((m-1)*3),1:3);
    rowIdx = bsxfun(@plus,((m-1)*NUM_FRAMES), 1:NUM_FRAMES);
    XYZ(rowIdx,:) = markerData.Pos(:,xyzIDX(idx));
end
xmin=floor(min(XYZ(:,1)))-500;  xmax=ceil(max(XYZ(:,1)))+500; 
ymin=floor(min(XYZ(:,2)))-500;  ymax=ceil(max(XYZ(:,2)))+500; 
zmin=floor(min(XYZ(:,3)))-500;  zmax=ceil(max(XYZ(:,3)))+500; 
% gmin = min([xmin ymin zmin]);
% gmax = max([xmax ymax zmax]);
% xlim([-2000 2000]);ylim([-2000 2000]);zlim([0 1500]);
xlim([xmin xmax]);ylim([ymin ymax]);zlim([zmin zmax]);
xlabel('Global x (mm)');ylabel('Global y (mm)');zlabel('Global z (mm)');
axis('square');view(az,el);grid on;
hTitle = title('Frame 1');

%% ---  initialise plot handles --- %
NUM_MARKERS = double(markerData.Names.Count);
hText       = zeros(NUM_MARKERS,1);
hMarkers    = zeros(NUM_MARKERS,1);
hLines      = [];
markerColor = [0.9 0.9 0];

%% ------------------------ Initialise Plot
% --- Plot marker names and dots --- %
hold on
keys = markerData.Names.keys;
for n=1:NUM_MARKERS
    markerName = keys{n};
    markerIdx  = markerData.Names(markerName);
    markerPos  = markerData.Pos(1,markerIdx);
    hText(n,1) = text(markerPos(1)+txtOffset,markerPos(2)+txtOffset,...
        markerPos(3)+txtOffset,markerName,'FontSize',8,'FontWeight','bold'); 
    hMarkers(n,1) = plot3(markerPos(1),markerPos(2),markerPos(3),...
        'Color',markerColor,'MarkerSize',2,'Marker','o','MarkerFaceColor',markerColor);
end
%% --- Plot lines joining dots together IF list provided
if bMarkerSet
    [NUM_PAIRS,~]=size(markerSet);
    hLines = zeros(NUM_PAIRS,1);
    for p=1:NUM_PAIRS
        pA = markerSet{p,1};
        pB = markerSet{p,2};
        posA = markerData.Pos(1,markerData.Names(pA));
        posB = markerData.Pos(1,markerData.Names(pB));
        if bMarkerColours
            hLines(p,1) = plot3([posA(1) posB(1)],[posA(2) posB(2)],...
            [posA(3) posB(3)],'LineWidth',2,'Color',markerSetColour{p});
        else
            hLines(p,1) = plot3([posA(1) posB(1)],[posA(2) posB(2)],...
                [posA(3) posB(3)],'LineWidth',2);
        end
    end
end
hold off;
% --- End plot --- %

%% -------------- package handles and data for updating figure
myData                = guihandles(hFig);
myData.hTitle         = hTitle;
myData.hLines         = hLines;
myData.hText          = hText;
myData.hMarkers       = hMarkers;
myData.txtOffset      = txtOffset; 
myData.markerData     = markerData;
myData.markerSet      = markerSet;
myData.bMarkerSet     = bMarkerSet;
myData.speedIncrement = speedIncrement;
myData.bSpin          = bSpin;
myData.NUM_FRAMES     = NUM_FRAMES;
myData.globalIdx      = globalIdx;
myData.playUntilFrame = globalIdx + DEFAULT_PLAYUNTILFRAMESKIP;

% initialise prev frame button
pushbuttonPrevHandle=uicontrol('Parent',hFig,'Style','PushButton',...
    'String','<<','Units','normalized','Position',[0.001 0.001 0.03 0.05]);
    set(pushbuttonPrevHandle,'Callback',@PushButtonPrevCallback);
myData.UI_CTRLS.pushbuttonPrevHandle=pushbuttonPrevHandle;

% initialise play/pause button
togbutPlayPauseHandle=uicontrol('Parent',hFig,'Style','togglebutton',...
    'String','|>','Units','normalized','Position',[0.031 0.001 0.03 0.05]);
    set(togbutPlayPauseHandle,'Callback',@PlayPauseCallback);
myData.UI_CTRLS.togbutPlayPauseHandle=togbutPlayPauseHandle;
    
% initialise next frame button
pushbuttonNextFrameHandle=uicontrol(...
    'Parent',hFig,'Style','PushButton',...
    'String','>>','Units','normalized','Position',[0.061 0.001 0.03 0.05]);
    set(pushbuttonNextFrameHandle,'Callback',@PushButtonNextCallback);
myData.UI_CTRLS.pushbuttonNextFrameHandle=pushbuttonNextFrameHandle;

% initialise editbox for video frame number
editVidFrameHandle = uicontrol('Parent',hFig,'Style','Edit',...
    'Units','normalized','Position',[0.091 0.001 0.075 0.05],...
    'String','Frame');
    set(editVidFrameHandle,'Callback',@EditFrameNumberCallBack);        
myData.UI_CTRLS.editVidFrameHandle=editVidFrameHandle;

editVidToFrameHandle = uicontrol('Parent',hFig,'Style','Edit',...
    'Units','normalized','Position',[0.181 0.001 0.075 0.05],...
    'String','Frame');
    set(editVidToFrameHandle,'Callback',@EditToFrameNumberCallBack);        
myData.UI_CTRLS.editVidToFrameHandle=editVidToFrameHandle;

guidata(hFig,myData);

end

function [myData] = updatePlotHandles(varargin)
if nargin == 2
    frameNum    = varargin{1};
    myData      = varargin{2};
end

markerData = myData.markerData;
markerSet  = myData.markerSet;
bMarkerSet = myData.bMarkerSet;
txtOffset  = myData.txtOffset; 
hText      = myData.hText;
hMarkers   = myData.hMarkers;
hTitle     = myData.hTitle;
hLines     = myData.hLines;

keys        = markerData.Names.keys;
NUM_MARKERS = markerData.Names.Count;

for n=1:NUM_MARKERS
    markerName = keys{n};
    markerIdx  = markerData.Names(markerName);
    markerPos  = markerData.Pos(frameNum,markerIdx);
    set(hText(n,1),'Position',[markerPos(1)+txtOffset,...
        markerPos(2)+txtOffset,markerPos(3)+txtOffset]); 
    set(hMarkers(n,1),'XDATA',markerPos(1),'YDATA',markerPos(2),...
        'ZDATA',markerPos(3));
end
% --- Plot lines joining dots together IF list provided
if bMarkerSet
    [NUM_PAIRS,~]=size(markerSet);
    for p=1:NUM_PAIRS
        pA = markerSet{p,1};
        pB = markerSet{p,2};
        posA = markerData.Pos(frameNum,markerData.Names(pA));
        posB = markerData.Pos(frameNum,markerData.Names(pB));
        set(hLines(p,1),'XDATA',[posA(1) posB(1)],...
            'YDATA',[posA(2) posB(2)],'ZDATA',[posA(3) posB(3)]);
    end
end    
set(hTitle,'String',['Frame ',num2str(frameNum)]);
[myData] = updateEditVidFrame(frameNum,myData);
%% ---- update data   
myData.markerData = markerData;
myData.markerSet  = markerSet;
myData.bMarkerSet = bMarkerSet;
myData.txtOffset  = txtOffset; 
myData.hText      = hText;
myData.hMarkers   = hMarkers;
myData.hTitle     = hTitle;
myData.hLines     = hLines;    
myData.globalIdx  = frameNum;
end

%PUSHBUTTONPREVCALLBACK callback function to handle loadvideo events
function PushButtonPrevCallback(hObject,event)
myData     = guidata(gcbo); % handle of callback whose object is executing
globalIdx  = myData.globalIdx;
prevIdx    = globalIdx-1;
myData     = checkValidFrame(prevIdx,myData);      
guidata(hObject,myData);
end

%PUSHBUTTONNEXTCALLBACK callback funcrtion to handle loadvideo events
function PushButtonNextCallback(hObject,event)
myData     = guidata(gcbo); % handle of callback whose object is executing
globalIdx  = myData.globalIdx;
nextIdx    = globalIdx +1;
[myData]   = checkValidFrame(nextIdx,myData);
guidata(hObject,myData);
end

%PLAYPAUSECALLBACK callback function to handle togglebutton events
% hObject - handles to object
function PlayPauseCallback(hObject,event)
myData = guidata(gcbo); % handle of callback whose object is executing
globalIdx  = myData.globalIdx;
NUM_FRAMES = myData.NUM_FRAMES;
    if get(hObject,'Value') == 0 % PAUSED
        set(hObject,'String','|>'); 
    else
        set(hObject,'String','||'); % PLAYING
        if isempty(globalIdx); 
            k=1;
        else
            k = globalIdx;
        end        
%         while k <NUM_FRAMES && get(hObject,'Value');
        while k < 1000 && get(hObject,'Value');
            [myData] = updatePlotHandles(k,myData);
            k=k+1;
            globalIdx=k;
            myData.globalIdx = globalIdx;
            guidata(hObject,myData);
            pause(0.0001);
        end         
        % if animation finished, need to stop
%         if get(hObject,'Value');
%             globalIdx = [];
%             myData.globalIdx = globalIdx;
%             guidata(hObject,myData);
%         end
    end        
    guidata(hObject,myData);
end

%EDITFRAMENUMBERCALLBACK callback function to handle userinput when the frame
%needs to be updated
function EditFrameNumberCallBack(hObject,event)
% checks that the input is an integer only.
% source: http://www.mathworks.com.au/matlabcentral/answers/14581-how-to-accept-only-numbers-in-a-edit-text-box
myData = guidata(gcbo); % handle of callback whose object is executing
str = get(hObject,'String');
frameNum = str2num(str);
[myData] = checkValidFrame(frameNum,myData);
guidata(hObject,myData);
end

%EDITTOFRAMENUMBERCALLBACK callback function to handle userinput when the frame
%needs to be updated
function EditToFrameNumberCallBack(hObject,event)
    % checks that the input is an integer only.
    % source: http://www.mathworks.com.au/matlabcentral/answers/14581-how-to-accept-only-numbers-in-a-edit-text-box
    myData = guidata(gcbo); % handle of callback whose object is executing
    str = get(hObject,'String');
    toFrameNum = str2num(str);
    if isempty(toFrameNum) % Check input is valid, i.e. number only
        set(hObject,'string','1');
        warndlg('Input must be numerical');
    elseif toFrameNum > myData.NUM_FRAMES || toFrameNum < 1
        warndlg(['Frame must be between, 1 & ',num2str(myData.NUM_FRAMES)]);
    else
        myData.playUntilFrame = toFrameNum;
    end
    guidata(hObject,myData);
end

% checks if the frame number is valid, if so updates the plot
function [myData]=checkValidFrame(frameNum,myData)
if isempty(frameNum) % Check input is valid, i.e. number only
    set(hObject,'string','1');
    warndlg('Input must be numerical');
elseif frameNum > myData.NUM_FRAMES || frameNum < 1
    warndlg(['Frame must be between, 1 & ',num2str(myData.NUM_FRAMES)]);
else
    myData = updatePlotHandles(frameNum,myData);
end
end

%UPDATEEDITFRAME updates the uicontrol edit object corresponding to the
%edit video frame box
function [myData] = updateEditVidFrame(frameNum,myData)
    set(myData.UI_CTRLS.editVidFrameHandle,'String',frameNum);
end




