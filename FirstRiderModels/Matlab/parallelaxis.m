function IB = parallelaxis(ICA, MA, a,b,c)
IB=ICA+MA*[b*b+c*c, -a*b, -a*c; -a*b, c*c+a*a, -b*c; -a*c, -b*c, a*a+b*b];