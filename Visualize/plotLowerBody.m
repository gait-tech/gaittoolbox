% ======================================================================
%> @brief Plot the lower body
%> @author Luke Sy
%>
%> Plot the lower body
%> 
%> @param pred predicted position
%> @param actl actual position
%> @param dim if 1 dimension equals (n x 3), else dimension equals (3 x n)
%>
%> @retval p plot object
% ======================================================================
function p = plotLowerBody(MP, LFEP, LFEO, LTIO,...
                           RFEP, RFEO, RTIO, marker)
    if nargin <= 7
        marker = '.';
    end
    
    validateattributes(MP, {'numeric'}, {});
    validateattributes(LFEP, {'numeric'}, {});
    validateattributes(LFEO, {'numeric'}, {});
    validateattributes(LTIO, {'numeric'}, {});
    validateattributes(RFEP, {'numeric'}, {});
    validateattributes(RFEO, {'numeric'}, {});
    validateattributes(RTIO, {'numeric'}, {});
    
    pelvL = line([LFEP(1) MP(1) RFEP(1)], [LFEP(2) MP(2) RFEP(2)],...
                 [LFEP(3) MP(3) RFEP(3)], ...
                 'Color', 'k', 'LineWidth', 2);
    pelvL.Marker = marker;
    lfemL = line([LFEP(1) LFEO(1)], [LFEP(2) LFEO(2)], ...
                 [LFEP(3) LFEO(3)],...
                 'Color', 'b', 'LineWidth', 2, 'LineStyle', '--');
    rfemL = line([RFEP(1) RFEO(1)], [RFEP(2) RFEO(2)], ...
                 [RFEP(3) RFEO(3)],...
                 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
    ltibL = line([LFEO(1) LTIO(1)], [LFEO(2) LTIO(2)], ...
                 [LFEO(3) LTIO(3)],...
                 'Color', 'b', 'LineWidth', 2);
    ltibL.Marker = marker;
    rtibL = line([RFEO(1) RTIO(1)], [RFEO(2) RTIO(2)], ...
                 [RFEO(3) RTIO(3)],...
                 'Color', 'r', 'LineWidth', 2);
    rtibL.Marker = marker;
end