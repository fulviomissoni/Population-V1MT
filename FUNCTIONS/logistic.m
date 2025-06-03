function CT = logistic(resp,logistic_centre,logistic_slope)

CT = 1./(1+exp(-logistic_slope*(resp-logistic_centre)));

end