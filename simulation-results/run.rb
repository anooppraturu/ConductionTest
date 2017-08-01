# run.rb: produces convergence data for anisotropic conduction
#
# usage: ruby run.rb
#
# from the analytic solution with a finite perpendicular conductivity,
# the central temperature reaches 1.01839 after 0.5 conduction times,
# defined as (sigma^2)/kappa_perp.
#
# the simulation dies when the central temperature reaches 1.01839,
# and we infer the effective perpendicular conduction time is twice
# the simulation time.
#

def issue_cmd(cmd)
  puts cmd
  system cmd
end

res   = [8, 16, 32, 64, 128]
tcond = Array.new

issue_cmd('mkdir -p out')
issue_cmd('mkdir -p err')
issue_cmd('mkdir -p numeric-profiles')

for i in res
  issue_cmd "./athena domain1/Nx1=#{i} domain1/Nx2=#{i} 1>out/#{i}.out 2>err/#{i}.err"
  issue_cmd "mv ic.dat  numeric-profiles/#{i}-ic.dat"
  issue_cmd "mv end.dat numeric-profiles/#{i}-end.dat"

  # simulation was designed to run for the equivalent of ~0.5
  # perpendicular conduction times... hence, we estimate the
  # conduction time as ~2 * the simulation run time
  #
  File.foreach("err/#{i}.err"){|x| tcond << 2.0 * x.split.last.to_f }
end


# print out convergence data
#
sigma     = 0.15
tcondpara = (sigma*sigma)/1.0

File.open('convergence.dat', 'w') do |file|
  printf(file,
         "# [1] = total res, [2] = tcond, [3] = (delta x)(grad ln T),  [4] = kappa_parallel/kappa_perp\n")
  res.each_index do |i|
    printf(file, "%6d\t%12.8f\t%6.3f\t%12.6f\n",
           res[i], tcond[i], res[i] * (sigma/2.0), tcond[i]/tcondpara)
  end
end

exit 0
