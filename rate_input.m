function axon = rate_input(rate, time_interval)
  axon = zeros(size(rate, 1), time_interval);
  one_num = round(rate * time_interval / 1000);
  for i = 1:size(rate, 1)
    if one_num(i) > time_interval
      k = 1:time_interval;
    else
      k = randperm(time_interval, one_num(i));
    end
    axon(i, k) = 1;
  end
end
