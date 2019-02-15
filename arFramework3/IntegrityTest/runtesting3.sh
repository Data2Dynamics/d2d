output=$(git pull)
if [ "$output" != "Already up-to-date." ]; then
  /usr/local/MATLAB/R2017a/bin/matlab -r "dotestserver3"&>matlab_output.txt
fi