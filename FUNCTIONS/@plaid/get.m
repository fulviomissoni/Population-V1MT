function value = get(p,type)

switch type
    case 'apert_rad'
        value = p.apert_rad;      
    case 'dur'
        value = p.dur;
    case 'truetheta'
        value = p.truetheta ;
    case 'vpld'
        value = p.vpld;
    case 'k'
        value = p.k;
    case 'vgrat'
        value = p.vgrat;
    case 'theta_g'
        value = p.theta_g;
    case 'alpha'
        value = p.alpha;
    case 'c'
        value = p.c;
    case 'mode'
        value = p.mode;
end