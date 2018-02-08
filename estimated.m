function xbar=estimated(input)

tau= input(1);
Ly= input(2:5);
in= input(6:9);
x= [in;tau];

xba= systemCTN(x);
xbar= xba+ Ly;
end
