function r = f(x, c)
  r = 0;
  powerOfX = 1;
  for i=1:length(c)
    r = r + c(i) * powerOfX;
    powerOfX = powerOfX * x;
  end