function beq = a2dR(d,n) 
% BEQ = A2DR(D, N) generates the decimal  
% equivalent beq of the binary representation  
% of a decimal number D with N bits for the  
% magnitude part obtained by rounding 
% 
m = 1; d1 = abs(d); 
while fix(d1) > 0 
	d1 = abs(d)/(10^m); 
	m = m+1; 
end 
beq = 0;d1 = d1 + 2^(-n-1); 
for k = 1:n 
	beq = fix(d1*2)/(2^k) + beq; 
	d1 = (d1*2) - fix(d1*2); 
end 
beq = sign(d).*beq.*10^(m-1); 