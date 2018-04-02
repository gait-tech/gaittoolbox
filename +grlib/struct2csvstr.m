% ======================================================================
%> @brief Converts struct to csv string
%>
%> @param s input struct
%> @param showcolname boolean if headers will be printed out (option)
%>
%> @retval out csv string
% ======================================================================
function out = struct2csvstr(s, showcolname)
    if nargin <= 1
        showcolname = false;
    end
    
    data = struct2table(s);
    out = grlib.table2csvstr(data, showcolname);
end