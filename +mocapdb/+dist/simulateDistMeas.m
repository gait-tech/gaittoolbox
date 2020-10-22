function out = simulateDistMeas(body, sigma, fname, israw)
	% Load CSV files
	%
    % :param body: input body
    % :type body: :class:`+pelib.@grBody`
    % :param sigma: gaussian noise level
	% :param fname: input filename
    % :type fname: String, Optional. If blank, won't output to file.
    % :param israw: does grBody only contain trimmed data
    % :type israw: Boolean, Optional. Defaults to false.
	%
	% :return: out - table containing distances
	%
	% .. Author: - Luke Sy (UNSW GSBME) 17 Jun 2020
    if nargin <= 1
        sigma = 0;
    end
    if nargin <= 2
        fname = '';
    end
    if nargin <= 3
        israw = false;
    end
    
    out = table;
    out.PV_LA = vecnorm(body.MIDPEL-body.LTIO, 2, 2) ...
        + normrnd(0, sigma, [body.nSamples, 1]);
    out.PV_RA = vecnorm(body.MIDPEL-body.RTIO, 2, 2) ...
        + normrnd(0, sigma, [body.nSamples, 1]);
    out.LA_RA = vecnorm(body.RTIO-body.LTIO, 2, 2) ...
        + normrnd(0, sigma, [body.nSamples, 1]);
    out.LLeg = vecnorm(body.LFEP-body.LTIO, 2, 2) ...
        + normrnd(0, sigma, [body.nSamples, 1]);
    out.RLeg = vecnorm(body.RFEP-body.RTIO, 2, 2) ...
        + normrnd(0, sigma, [body.nSamples, 1]);

    if ~isempty(fname)
        fid = fopen(fname, 'w');
        fprintf(fid, "// IsRaw=%d Version=%.2f\n", israw, 1.0);
        fclose(fid);
        writetable(out, fname, 'WriteMode', 'append', ...
                   'WriteVariableNames', true);
    end
end

