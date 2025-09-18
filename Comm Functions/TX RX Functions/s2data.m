function data = s2data(s_mod,S_alphabet)

sizes = size(s_mod);
[~,orient] = max(sizes);

cmp_matrix = abs(s_mod.' - S_alphabet).^2;

[~,data] = min(cmp_matrix,[],orient);
data = data.';


end