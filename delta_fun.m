function delta = delta_fun(x)
  delta = zeros(1, length(x));
  i = find(x == 0);
  delta(i) = 1;
end
