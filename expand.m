% Expand column vector v to d columns
function m = expand(v, d)
  m = zeros(length(v), d);
  for i = 1:d
    m(:, i) = v;
  end
end
