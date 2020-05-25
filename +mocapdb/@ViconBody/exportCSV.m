function exportCSV(obj, fname, info)
	% Export as CSV file
    % 
	% :param obj: this ViconBody
	% :param fname: filename of file to be saved
    % :param info: additional info to be added at each csv file
	%
	% .. Author: - Luke Sy (UNSW GSBME) 2020 May 23

    if nargin <= 2
        info = "";
    end
    
    info0 = sprintf("// fs=%d  nSamples=%d\n// frame=%s\n// posUnit=%s\n// origFName=%s", ...
                        obj.fs, obj.nSamples, obj.frame, obj.posUnit, obj.srcFileName);
    
    %% store object attributes to struct
    out = struct();
    for i=obj.posList
        out.(i{1}) = obj.(i{1});
    end
    for i=obj.oriList
        out.(i{1}) = obj.(i{1});
    end
    out = struct2table(out);
    
    %% save struct and info to csv
    if ~endsWith(fname, '.csv')
        fnameFinal = sprintf("%s.csv", fname);
    end
    fid = fopen(fnameFinal, 'w');
    if fid < 0
        error("Can't open %s", fnameFinal);
    end
    fprintf(fid, "%s\n", info0);
    fprintf(fid, "// %s\n", info);
    fclose(fid);
    writetable(out, fnameFinal, 'WriteMode', 'append', ...
               'WriteVariableNames', true);
end