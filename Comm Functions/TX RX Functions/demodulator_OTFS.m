function data_demod = demodulator_OTFS(data, M)
% function data_demod = demod_mqam(data, M)
%    M: 4, 16 or 64
%

Es = 1; % the energy per symbol
% the size of one decision area 2*delta x 2*delta
delta = sqrt(3*Es/(2*(M-1)));

bignumber = 1000;

N_data = length(data);

if M == 2
   data_demod = zeros(1, length(data));
   %data_demod(real(data) < 0) = 1;
   data_demod(real(data) > 0) = 1;
elseif M == 4
   data_demod = zeros(1, 2*N_data);
   data_demod_inphase = zeros(1, N_data);
   data_demod_quadrature = zeros(1, N_data);

   amplitude_level = 0;
   
   data_sym_inphase = real(data);
   data_sym_quadrature = imag(data);
   mapping = [0 1];  

   level_idx_inphase = (data_sym_inphase < 0);
   data_demod_inphase(level_idx_inphase) = mapping(2);
   
   level_idx_quadrature = (data_sym_quadrature < 0);
   data_demod_quadrature(level_idx_quadrature) = mapping(2);

   data_demod(1:2:end) = data_demod_inphase;
   data_demod(2:2:end) = data_demod_quadrature;
elseif M == 8
   data_demod = demod_8psk(data);
elseif M == 16
   data_demod = zeros(1, 4*N_data);
   data_demod_inphase_1 = zeros(1, N_data);
   data_demod_inphase_2 = zeros(1, N_data);
   data_demod_quadrature_1 = zeros(1, N_data);
   data_demod_quadrature_2 = zeros(1, N_data);

   amplitude_level = [-2:2:2]*delta;
   N_level = length(amplitude_level);
   
   data_sym_inphase = real(data);
   data_sym_quadrature = imag(data);

   mapping_1 = [1 1 0 0];
   mapping_2 = [1 0 0 1];
   flag_inphase = logical(zeros(1, N_data));
   flag_quadrature = logical(zeros(1, N_data));
   for k = 1:N_level
      level_idx_inphase = (data_sym_inphase <= amplitude_level(k));
      data_sym_inphase(level_idx_inphase) = bignumber;
      flag_inphase(level_idx_inphase) = logical(1);
      data_demod_inphase_1(level_idx_inphase) = mapping_1(k);
      data_demod_inphase_2(level_idx_inphase) = mapping_2(k);
      
      level_idx_quadrature = (data_sym_quadrature <= amplitude_level(k));
      data_sym_quadrature(level_idx_quadrature) = bignumber;
      flag_quadrature(level_idx_quadrature) = logical(1);
      data_demod_quadrature_1(level_idx_quadrature) = mapping_1(k);
      data_demod_quadrature_2(level_idx_quadrature) = mapping_2(k);      
   end   
   
   data_demod_inphase_1(~flag_inphase) = mapping_1(N_level+1);
   data_demod_inphase_2(~flag_inphase) = mapping_2(N_level+1);
   data_demod_quadrature_1(~flag_quadrature) = mapping_1(N_level+1);
   data_demod_quadrature_2(~flag_quadrature) = mapping_2(N_level+1);
   
   
   data_demod(1:4:end) = data_demod_inphase_1;
   data_demod(3:4:end) = data_demod_inphase_2;
   data_demod(2:4:end) = data_demod_quadrature_1;
   data_demod(4:4:end) = data_demod_quadrature_2;

elseif M == 64

   data_demod = zeros(1, 6*N_data);
   data_demod_inphase_1 = zeros(1, N_data);
   data_demod_inphase_2 = zeros(1, N_data);
   data_demod_inphase_3 = zeros(1, N_data);
   data_demod_quadrature_1 = zeros(1, N_data);
   data_demod_quadrature_2 = zeros(1, N_data);
   data_demod_quadrature_3 = zeros(1, N_data);

   amplitude_level = [-6:2:6]*delta;
   N_level = length(amplitude_level);
   
   data_sym_inphase = real(data);
   data_sym_quadrature = imag(data);

   mapping_1 = [1 1 1 1 0 0 0 0];
   mapping_2 = [1 1 0 0 0 0 1 1];
   mapping_3 = [1 0 0 1 1 0 0 1];
   flag_inphase = logical(zeros(1, N_data));
   flag_quadrature = logical(zeros(1, N_data));
   for k = 1:N_level
      level_idx_inphase = (data_sym_inphase <= amplitude_level(k));
      data_sym_inphase(level_idx_inphase) = bignumber;
      flag_inphase(level_idx_inphase) = logical(1);
      data_demod_inphase_1(level_idx_inphase) = mapping_1(k);
      data_demod_inphase_2(level_idx_inphase) = mapping_2(k);
      data_demod_inphase_3(level_idx_inphase) = mapping_3(k);
      
      level_idx_quadrature = (data_sym_quadrature <= amplitude_level(k));
      data_sym_quadrature(level_idx_quadrature) = bignumber;
      flag_quadrature(level_idx_quadrature) = logical(1);
      data_demod_quadrature_1(level_idx_quadrature) = mapping_1(k);
      data_demod_quadrature_2(level_idx_quadrature) = mapping_2(k);      
      data_demod_quadrature_3(level_idx_quadrature) = mapping_3(k);      
   end   
   
   data_demod_inphase_1(~flag_inphase) = mapping_1(N_level+1);
   data_demod_inphase_2(~flag_inphase) = mapping_2(N_level+1);
   data_demod_inphase_3(~flag_inphase) = mapping_3(N_level+1);
   data_demod_quadrature_1(~flag_quadrature) = mapping_1(N_level+1);
   data_demod_quadrature_2(~flag_quadrature) = mapping_2(N_level+1);
   data_demod_quadrature_3(~flag_quadrature) = mapping_3(N_level+1);
      
   data_demod(1:6:end) = data_demod_inphase_1;
   data_demod(3:6:end) = data_demod_inphase_2;
   data_demod(5:6:end) = data_demod_inphase_3;
   data_demod(2:6:end) = data_demod_quadrature_1;
   data_demod(4:6:end) = data_demod_quadrature_2;
   data_demod(6:6:end) = data_demod_quadrature_3;   
elseif M == 256 
   data_demod = zeros(1, 8*N_data);
   data_demod_inphase_1 = zeros(1, N_data);
   data_demod_inphase_2 = zeros(1, N_data);
   data_demod_inphase_3 = zeros(1, N_data);
   data_demod_inphase_4 = zeros(1, N_data);
   data_demod_quadrature_1 = zeros(1, N_data);
   data_demod_quadrature_2 = zeros(1, N_data);
   data_demod_quadrature_3 = zeros(1, N_data);
   data_demod_quadrature_4 = zeros(1, N_data);

   amplitude_level = [-14:2:14]*delta;
   N_level = length(amplitude_level);
   
   data_sym_inphase = real(data);
   data_sym_quadrature = imag(data);

   mapping_1 = [1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0];
   mapping_2 = [0 0 0 0 1 1 1 1 1 1 1 1 0 0 0 0];
   mapping_3 = [1 1 0 0 0 0 1 1 1 1 0 0 0 0 1 1];
   mapping_4 = [1 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1];
   flag_inphase = logical(zeros(1, N_data));
   flag_quadrature = logical(zeros(1, N_data));
   for k = 1:N_level
      level_idx_inphase = (data_sym_inphase <= amplitude_level(k));
      data_sym_inphase(level_idx_inphase) = bignumber;
      flag_inphase(level_idx_inphase) = logical(1);
      data_demod_inphase_1(level_idx_inphase) = mapping_1(k);
      data_demod_inphase_2(level_idx_inphase) = mapping_2(k);
      data_demod_inphase_3(level_idx_inphase) = mapping_3(k);
      data_demod_inphase_4(level_idx_inphase) = mapping_4(k);
      
      level_idx_quadrature = (data_sym_quadrature <= amplitude_level(k));
      data_sym_quadrature(level_idx_quadrature) = bignumber;
      flag_quadrature(level_idx_quadrature) = logical(1);
      data_demod_quadrature_1(level_idx_quadrature) = mapping_1(k);
      data_demod_quadrature_2(level_idx_quadrature) = mapping_2(k);      
      data_demod_quadrature_3(level_idx_quadrature) = mapping_3(k);      
      data_demod_quadrature_4(level_idx_quadrature) = mapping_4(k);      
   end   
   
   data_demod_inphase_1(~flag_inphase) = mapping_1(N_level+1);
   data_demod_inphase_2(~flag_inphase) = mapping_2(N_level+1);
   data_demod_inphase_3(~flag_inphase) = mapping_3(N_level+1);
   data_demod_inphase_4(~flag_inphase) = mapping_4(N_level+1);
   data_demod_quadrature_1(~flag_quadrature) = mapping_1(N_level+1);
   data_demod_quadrature_2(~flag_quadrature) = mapping_2(N_level+1);
   data_demod_quadrature_3(~flag_quadrature) = mapping_3(N_level+1);
   data_demod_quadrature_4(~flag_quadrature) = mapping_4(N_level+1);
      
   data_demod(1:8:end) = data_demod_inphase_1;
   data_demod(3:8:end) = data_demod_inphase_2;
   data_demod(5:8:end) = data_demod_inphase_3;
   data_demod(7:8:end) = data_demod_inphase_4;
   data_demod(2:8:end) = data_demod_quadrature_1;
   data_demod(4:8:end) = data_demod_quadrature_2;
   data_demod(6:8:end) = data_demod_quadrature_3;   
   data_demod(8:8:end) = data_demod_quadrature_4;   
end

  
