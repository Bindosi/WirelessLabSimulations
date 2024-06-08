x = ones(1e2,1)

s.T ="RRC"
s.sps = 10
s.span = 10
s.alpha = 0.4

fltr(x, s)