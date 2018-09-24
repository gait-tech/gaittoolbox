qSenseLTIB = qOri.sv.LTIB;
qViconLTIB = qOri.v.LTIB;
qSenseLTIB2 = qalign(qSenseLTIB, qViconLTIB);

qSenseRTIB = qOri.sv.RTIB;
qViconRTIB = qOri.v.RTIB;
qSenseRTIB2 = qalign(qSenseRTIB, qViconRTIB);

clf; pelib.viz.plotQuaternion(qSenseLTIB2, qViconLTIB);
pelib.viz.plotQuaternion(qSenseRTIB2, qViconRTIB);

function qNew = qalign(q, qRef)
    flip = xor(q(:,1) > 0, qRef(:,1) > 0);
    qNew = q;
    qNew(flip, :) = q(flip, :) .* -1;
end