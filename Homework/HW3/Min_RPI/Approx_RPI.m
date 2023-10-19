function F_approx = Approx_RPI(A,W,epsilon)

if isprop(W,'H')
    % Choose initial s
    s = 0;
    % Set maximum iterations
    for i = 1:100
        % Increment s
        s = s + 1
        % Compute alpha min of s
        alpha_o = alpha_min(W,A,s);
        % Compute M(s)
        M_s = M_of_s(W,A,s);
        if alpha_o <= epsilon/(epsilon+M_s);
            break
        end
    end
    
    Fs = W;
    for i = 1:s % because nilpotent in one step
        Fs = plus(Fs,affineMap(W,A^i));
        Fs.minHRep;
    end
    F_approx = 1/(1-alpha_o)*Fs;
    F_approx.minHRep;
    % F_approx.plot
elseif isfield(W,'G')
    % Choose initial s
    s = 0;
    % Set maximum iterations
    for i = 1:100
        % Increment s
        s = s + 1
        % Compute alpha min of s
        alpha_o = alpha_min(W,A,s);
        % Compute M(s)
        M_s = M_of_s(W,A,s);
        if alpha_o <= epsilon/(epsilon+M_s);
            break
        end
    end
    
    Fs = W;
    for i = 1:s % because nilpotent in one step
        Fs.c = Fs.c + A^i*W.c;
        Fs.G = [Fs.G A^i*W.G];
    end
    F_approx.c = 1/(1-alpha_o)*Fs.c;
    F_approx.G = 1/(1-alpha_o)*Fs.G;
else
    disp('Error in Approx_RPI.m: Set must be in either H-Rep or G-Rep')
    F_approx = [];
end