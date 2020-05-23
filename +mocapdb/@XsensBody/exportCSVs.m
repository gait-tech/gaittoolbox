function exportCSVs(obj, fname, info)
	% Export all sensor information as CSV files
    % 
    % File name convension: fname-bodypart.csv
    % Each file will contain acceleration, gyroscope, magnetometer, and 
    % quaternion data
	%
	% :param obj: this XsensBody or struct('bodypart', table)
	% :param fname: filename of file to be saved
    % :param info: additional info to be added at each csv file
	%
	% .. Author: - Luke Sy (UNSW GSBME) 2020 May 22

    if nargin <= 2
        info = "";
    end
    
    info0 = sprintf("// fs=%d nSamples=%d\n// frame=%s\n// origFName=%s", ...
                        obj.fs, obj.nSamples, obj.frame, obj.srcFileName);
    
    for i=obj.segList
        if isprop(obj, i{1}) && ~isempty(obj.(i{1}))
            fnameFinal = sprintf("%s-%s", fname, i{1});
            if ~endsWith(fnameFinal, '.csv')
                fnameFinal = sprintf("%s.csv", fnameFinal);
            end
            fid = fopen(fnameFinal, 'w');
            if fid < 0
                error("Can't open %s", fnameFinal);
            end
            fprintf(fid, "%s\n", info0);
            fprintf(fid, "// %s\n", info);
            fclose(fid);
            writetable(obj.(i{1}), fnameFinal, 'WriteMode', 'append', ...
                       'WriteVariableNames', true);
        end
    end
end