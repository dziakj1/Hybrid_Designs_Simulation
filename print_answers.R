# Print answers in format of paper;
load("did-simulations.rdata"); 
print(round(t(proximal_T1E),4));
print(round(t(distal_T1E),4));
print(round(t(proximal_power),4));
print(round(t(distal_power),4));