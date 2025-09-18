function RzDD = make_RzDD(obj)

% % Debugging inputs
% Fd = 15e3;
% shape = "sinc";
% N_tsyms = 2;
% M_subcarriers = 8;
% N0 = 1;

% Inputs from obj
subcarrier_spacing = 1 / obj.T1;
select_filter = obj.select_filter;
N_tsyms = obj.N_tsyms;
M_sbcars = obj.M_sbcars;
N0 = obj.N0;

% Internal settings
T = 1/subcarrier_spacing;

% Create each element of RzTF
RzTF = zeros(N_tsyms*M_sbcars);
for m1 = 0:M_sbcars-1
    for m2 = 0:M_sbcars-1
        for n1 = 0:N_tsyms-1
            for n2 = 0:N_tsyms-1
                RzTF(m1*N_tsyms+n1+1,m2*N_tsyms+n2+1) = N0 * ...
                    exp(1j*2*pi*m2*(n1-n2)*subcarrier_spacing*T) * ...
                    xambig((n1-n2)*T,(m1-m2)*subcarrier_spacing,T,select_filter);
            end
        end
    end
end

% Create DFT matrices and needed Kronecker product
F_N = fr_DFT(N_tsyms);
F_M = fr_DFT(M_sbcars);
F_cur = kron(F_M,F_N');

% Create final result
RzDD = F_cur' * RzTF * F_cur;