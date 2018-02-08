function xbar=estimated_rand(input)
tau= input(1);
Ly= input(2:5);
in= input(6:9);
x= [in;tau];

xba= systemCTN_rand(x);
xbar= xba+ Ly;
end
