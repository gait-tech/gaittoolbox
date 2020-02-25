function [ hfig ] = updateFigureContents( figstr )
%UPDATEFIGURECONTENTS Retrieves the handle of an existing figure with name
%'figstr'
%   Detailed explanation goes here
[fignames,hfigs] = getFigNames;
fighandle        = findFigureHandle(fignames,figstr );
if isempty(fighandle)
    hfig = figure('name',figstr,'Color','w');
else
    hfig = figure(hfigs(fighandle));clf
end

end

function [ fnames,hfigs ] = getFigNames(  )
%GETFIGNAMES Summary of this function goes here
%   Detailed explanation goes here
hfigs  = findobj('Type','figure');
fnames = get(hfigs,'Name');
if (length(hfigs) == 1)
    fnames = {fnames};
end

end

function [ fighandle ] = findFigureHandle(fignames,figname )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

fighandle = find(strcmp(fignames,figname)==1);
            
end