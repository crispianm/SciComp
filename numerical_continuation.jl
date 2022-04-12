function natural_parameter_continuation(f, u0_, alpha_range, delta_alpha, discretisation, arg...)
    
    num_alphas = floor(abs(alpha_range[end] - alpha_range[1]) / delta_alpha)
    alphas = linrange(alpha_range[1], alpha_range[end], num_alphas)
    return alphas

end

function pseudo_arclength(f, u0_, alpha_range, delta_alpha, discretisation, arg...)

    

end
