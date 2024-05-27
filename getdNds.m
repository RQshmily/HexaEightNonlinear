function [ dnds ] = getdNds(ksi,eta,zeta)

dnds=0.125*[
    -(1-eta)*(1-zeta),(1-eta)*(1-zeta),(1+eta)*(1-zeta),-(1+eta)*(1-zeta),...
    -(1-eta)*(1+zeta),(1-eta)*(1+zeta),(1+eta)*(1+zeta),-(1+eta)*(1+zeta);

    -(1-ksi)*(1-zeta),-(1+ksi)*(1-zeta),(1+ksi)*(1-zeta),(1-ksi)*(1-zeta),...
    -(1-ksi)*(1+zeta),-(1+ksi)*(1+zeta),(1+ksi)*(1+zeta),(1-ksi)*(1+zeta);

    -(1-ksi)*(1-eta),-(1+ksi)*(1-eta),-(1+ksi)*(1+eta),-(1-ksi)*(1+eta),...
    (1-ksi)*(1-eta), (1+ksi)*(1-eta), (1+ksi)*(1+eta), (1-ksi)*(1+eta);
    ];

end