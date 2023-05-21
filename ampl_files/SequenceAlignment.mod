
param n;                               # length of sequence s1
param m;                               # length of sequence s2
param mis_p;                           # mismatch penality
param gap_p;                           # gap penality
param r;                               # match reward
param a{i in 1..n, j in 1..m} binary;  # matrix that indicates whether s1[i] is equal to s2[j]

var x{i in 1..n, j in 1..m} binary;    # matrix that indicates whether s1[i] and s2[j] are aligned
var y{i in 1..n} binary;               # array that indicates whether s1[i] is aligned with a gap in s2
var z{j in 1..m} binary;               # array that indicates whether s2[j] is aligned with a gap in s1
var v;

# a symbol in a sequence can't be aligned with both another symbol and with a gap at the same time.
subject to maxAlignmentsS1{i in 1..n}: 
	(sum{j in 1..m} x[i,j]) + y[i] <= 1;
subject to maxAlignmentsS2{j in 1..m}: 
	(sum{i in 1..n} x[i,j]) + z[j] <= 1;
    
# if s1[i] is aligned with s2[j], previous symbols in s1 can't be aligned with following symbols 
# in s2, and vice versa.
subject to alignmentOrder1{i in 1..n, j in 1..m, k in 1..i-1, h in j+1..m}:
    x[i,j] + x[k,h] <= 1;
subject to alignmentOrder2{i in 1..n, j in 1..m, k in i+1..n, h in 1..j-1}:
    x[i,j] + x[k,h] <= 1;
    
# two symbols to be aligned must be the same.
subject to sameSymbol{i in 1..n, j in 1..m}: 
	x[i,j] <= a[i,j];
    
# if two symbols are aligned but have different positions in their sequence, there must be 
# enough gaps in previous positions of both the sequences.
subject to previousGaps{i in 1..n, j in 1..m}:
	x[i,j] * ( ( sum{k in 1..i} y[k] ) - ( sum{h in 1..j} z[h] ) ) == x[i,j] * (i - j);

# to know the number of final gaps in the shortest sequence, we need to calculate the absolute value 
# of the difference between the two aligned sequences. Since the absolute value is not a linear function,
# we add the following two constraints to set the variable v to be the absolute value we need.
subject to absoluteValueA:
	-v <= (n + sum{j in 1..m} z[j]) - (m + sum{i in 1..n} y[i]);
subject to absoluteValueB:
	(n + sum{j in 1..m} z[j]) - (m + sum{i in 1..n} y[i]) <= v;
    
# with this objective function we try to maximize the number of matches and 
# at the same time minimize the number of gaps and mismatches.
minimize obj: 
	( sum{j in 1..m} sum{i in 1..n} (x[i,j]) * r) + #matches
	( sum{i in 1..n} y[i] * gap_p )               + #gaps
	( sum{j in 1..m} z[j] * gap_p )               + 
	( sum{i in 1..n} (1 - ( ( sum{j in 1..m} x[i,j] ) + y[i] ) ) * mis_p ) +  #mismatches
	( v * gap_p ); #ending gaps
	
