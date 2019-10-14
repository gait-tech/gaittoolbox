function out = diff(obj1, obj2, seq)
	% Returns the difference between grBody1 and grBody2
	%
	% :param obj1: grBody 1 (self)
	% :param obj2: grBody 2 (other)
	% :param seq: orientation sequence
	%
	% :return: out - struct with the difference of pos and ori parameters
	%
	% .. Author: - Luke Sy (UNSW GSBME)

    if nargin <= 2
        seq = 'YXZ';
    end
    
    out = struct;
    
    posList = obj1.posList;
    for i=1:length(posList)
        if length(obj1.(posList{i})) == 0 || length(obj2.(posList{i})) == 0
            out.(posList{i}) = [];
        else
            out.(posList{i}) = obj1.(posList{i}) - obj2.(posList{i});
        end
    end
    
    [r1 r2 r3] = quat2angle(obj1.qRPV, seq);
    eul1 = [r1 r2 r3];
    [r1 r2 r3] = quat2angle(obj2.qRPV, seq);
    eul2 = [r1 r2 r3];
    out.qRPV = (eul1 - eul2)*180/pi;   
    out.qLHIP = (obj1.calcJointAnglesLHip() - obj2.calcJointAnglesLHip())*180/pi;
    out.qRHIP = (obj1.calcJointAnglesRHip() - obj2.calcJointAnglesRHip())*180/pi;
    out.qLKNE = (obj1.calcJointAnglesLKnee() - obj2.calcJointAnglesLKnee())*180/pi;
    out.qRKNE = (obj1.calcJointAnglesRKnee() - obj2.calcJointAnglesRKnee())*180/pi;
    
%     oriList = obj1.oriList;
%     for i=1:length(oriList)
%         if length(obj1.(oriList{i})) == 0 || length(obj2.(oriList{i})) == 0
%             out.(oriList{i}) = [];
%         else
%             [r1 r2 r3] = quat2angle(obj1.(oriList{i}), seq);
%             eul1 = [r1 r2 r3];
%             [r1 r2 r3] = quat2angle(obj2.(oriList{i}), seq);
%             eul2 = [r1 r2 r3];
% 
%             if obj1.oriUnit == 'deg'
%                 out.(oriList{i}) = pelib.anglediff(eul1, eul2)*180/pi;
%             else
%                 out.(oriList{i}) = pelib.anglediff(eul1, eul2);
%             end
%         end
%     end
end