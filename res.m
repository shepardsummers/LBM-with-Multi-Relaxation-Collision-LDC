function [out, max_error] = res(A, B, min_error)
    
    error = A - B;
    max_error = max(error, [], "all");

    %disp(m)ax_error

    if max_error < min_error
        out = false;
    else
        out = true;
    end

end

