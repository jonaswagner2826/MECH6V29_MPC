function M_s = M_of_s(W,A,s)

if isprop(W,'H')
    I = eye(size(A,1));
    sum_max = zeros(length(I),1);
    for j = 1:length(I)
        ej = I(:,j);
        sum_plus = 0;
        sum_minus = 0;
        for i = 0:s-1
            a = (A^i)'*ej;
            hW_a_plus = support(W,a);
            hW_a_minus = support(W,-a);
            sum_plus = sum_plus + hW_a_plus;
            sum_minus = sum_minus + hW_a_minus;
        end
        sum_max(j) = max(sum_plus,sum_minus);
    end
    M_s = max(sum_max);
elseif isfield(W,'G')
    I = eye(size(A,1));
    sum_max = zeros(length(I),1);
    for j = 1:length(I)
        ej = I(:,j);
        sum_plus = 0;
        sum_minus = 0;
        for i = 0:s-1
            a = (A^i)'*ej;
            hW_a_plus = sum(abs(a'*W.G)) + a'*W.c;
            hW_a_minus = sum(abs(a'*W.G)) - a'*W.c;
            sum_plus = sum_plus + hW_a_plus;
            sum_minus = sum_minus + hW_a_minus;
        end
        sum_max(j) = max(sum_plus,sum_minus);
    end
    M_s = max(sum_max);
else
    disp('Error in M_of_s.m: Set must be in either H-Rep or G-Rep')
    M_s = [];
end
