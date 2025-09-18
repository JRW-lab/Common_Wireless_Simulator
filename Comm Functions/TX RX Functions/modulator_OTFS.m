function [data_mod, group0, group1] = modulator_OTFS(data, M, group0, group1)
% function data_mod = modulator(data, M, group0, group1)
%    M: 2, 4, 8, 16 or 64
%    the output of the demodulator are soft values, group0 and group1 are
%    used by the demodulator to generate these soft values
%

if nargin == 2
   group0 = 1;
   group1 = 1;
end


% convert -1 to 1
data(data<0) = 0;

% used to indicate if the modulator needs to generate group0 and group1
group_flag = sum(group0);

Es = 1; % the energy per symbol
% the size of one decision area 2*delta x 2*delta
delta = sqrt(3*Es/(2*(M-1)));

% N_data = length(data);

% bpsk is not part of HSDPA, I add it for compatible issue
if M == 2
   %mapping = [1, -1];
   mapping = [-1, 1];
   data_mod = mapping(data+1);
elseif M == 4
   mapping = [1, -1]*delta;

   data_inphase = data(1:2:end);
   data_sym_inphase = mapping(data_inphase+1);

   data_quadrature = data(2:2:end);
   data_sym_quadrature = mapping(data_quadrature+1);

   data_mod = data_sym_inphase+1j*data_sym_quadrature;
elseif M == 8
% if 0
%    mapping_idx = [0 1 3 2 7 6 4 5];
%    mapping = exp(j*(mapping_idx*pi/4+pi/8));
%    data_1 = data(1:3:end);
%    data_2 = data(2:3:end);
%    data_3 = data(3:3:end);
%    data_mod = mapping(data_3+data_2*2+data_1*4+1);
% end
   if ~isempty(data)
      data_mod = eightpsk(data);
   else
      data_mod = zeros(0,0);
   end
elseif M == 16
   mapping = [1, 3, -1, -3]*delta;

   data_inphase_1 = data(1:4:end);
   data_inphase_2 = data(3:4:end);
   data_sym_inphase = mapping((data_inphase_2 + data_inphase_1*2)+1);
   
   data_quadrature_1 = data(2:4:end);
   data_quadrature_2 = data(4:4:end);
   data_sym_quadrature = mapping((data_quadrature_2 + data_quadrature_1*2)+1);

   data_mod = data_sym_inphase+1j*data_sym_quadrature;
elseif M == 64
   mapping = [3 1 5 7 -3 -1 -5 -7]*delta;
   data_inphase_1 = data(1:6:end);
   data_inphase_2 = data(3:6:end);
   data_inphase_3 = data(5:6:end);
   data_sym_inphase = mapping((data_inphase_3+data_inphase_2*2+data_inphase_1*4)+1);
  
   data_quadrature_1 = data(2:6:end);
   data_quadrature_2 = data(4:6:end);
   data_quadrature_3 = data(6:6:end);
   data_sym_quadrature = mapping((data_quadrature_3+data_quadrature_2*2+data_quadrature_1*4)+1);

   data_mod = data_sym_inphase+1j*data_sym_quadrature;
elseif M == 256
   mapping = [11 9 13 15 5 7 3 1 -11 -9 -13 -15 -5 -7 -3 -1]*delta;
   data_inphase_1 = data(1:8:end);
   data_inphase_2 = data(3:8:end);
   data_inphase_3 = data(5:8:end);
   data_inphase_4 = data(7:8:end);
   data_sym_inphase = mapping((data_inphase_4+data_inphase_3*2+data_inphase_2*4+data_inphase_1*8)+1);
  
   data_quadrature_1 = data(2:8:end);
   data_quadrature_2 = data(4:8:end);
   data_quadrature_3 = data(6:8:end);
   data_quadrature_4 = data(8:8:end);
   data_sym_quadrature = mapping((data_quadrature_4+data_quadrature_3*2+data_quadrature_2*4+data_quadrature_1*8)+1);

   data_mod = data_sym_inphase+1j*data_sym_quadrature;
end

if group_flag == 0
  count_0 = ones(1, log2(M));
  count_1 = ones(1, log2(M));
  for k = 1:M
      k_bin = dec2bin(k-1, log2(M))-'0';
      % recursive call itself, to generate the modualtion symbol of k-1
      k_sym = modulator_OTFS(k_bin, M, 1, 1); 
      for t = 1:log2(M)
         if (k_bin(t) == 0) % the t-th bit of k_bin is 0
            group0(t, count_0(t)) = k_sym;
            count_0(t) = count_0(t)+1;
         else   % the t-th bit of k_bin is 1
            group1(t, count_1(t)) = k_sym;
            count_1(t) = count_1(t)+1;
         end
      end
   end     
end

