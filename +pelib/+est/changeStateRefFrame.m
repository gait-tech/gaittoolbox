% .. Author: - Luke Sy (UNSW GSBME) - 9/22/18
% @brief change reference frame of state
% Supported changes are vicon -> MIDPEL
%
% :param state input kf_3_kmus_v3 state
% :param ref reference frame
% :param format state format (kf_3_kmus_v3)
%
% :return: out state at new frame
% .. Author: - Luke Sy (UNSW GSBME) - 9/22/18
function out = changeStateRefFrame(state, ref, format)
    if nargin <= 1
        ref = 'MIDPEL';
        format = 'kf_3_kmus_v3';
    end
           
    if strcmp(format, 'kf_3_kmus_v3')
        refPos = 1:3;
        refVel = 4:6;
        refOri = 7:10;
        posList = [1:3 11:13 21:23];
        velList = [4:6 14:16 24:26];
        oriList = [7:10 17:20 27:30];
    end
    
    if isstruct(state)
        out = state;
        sList = {'predState', 'zuptState', 'cstrState'};
        for j=1:length(sList)
            sName = sList{j};
            for i=1:3:length(posList)
                target = posList(i:i+2);
                out.(sName)(:,target) = quatrotate(state.(sName)(:,refOri), ...
                    state.(sName)(:,target)-state.(sName)(:,refPos));
            end

            for i=1:3:length(velList)
                target = velList(i:i+2);
                out.(sName)(:,target) = quatrotate(state.(sName)(:,refOri), ...
                    state.(sName)(:,target)-state.(sName)(:,refVel));
            end

            for i=1:4:length(oriList)
                target = oriList(i:i+3);
                out.(sName)(:,target) = quatmultiply(quatconj(state.(sName)(:,refOri)), ...
                    state.(sName)(:,target));
            end
        end
    else
        out = state;
        for i=1:3:length(posList)
            target = posList(i:i+2);
            out(:,target) = quatrotate(state(:,refOri), ...
                state(:,target)-state(:,refPos));
        end

        for i=1:3:length(velList)
            target = velList(i:i+2);
            out(:,target) = quatrotate(state(:,refOri), ...
                state(:,target)-state(:,refVel));
        end

        for i=1:4:length(oriList)
            target = oriList(i:i+3);
            out(:,target) = quatmultiply(quatconj(state(:,refOri)), ...
                state(:,target));
        end
    end
end