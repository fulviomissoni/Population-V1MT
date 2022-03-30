function set(p,type,value)

%eval(['p.',type,'=',value,';']);

switch type
        case 'apert_rad'
        p.apert_rad = value;      
    case 'dur'
        p.dur = value;
    case 'truetheta'
        p.truetheta = value;
    case 'vpld'
        p.vpld = value;
    case 'k'
        p.k = value;
    case 'vgrat'
        p.vgrat = value;
    case 'theta_g'
        p.theta_g = value;
    case 'alpha'
        p.alpha = value;
    case 'c'
        p.c = value;
    case 'mode'
        p.mode = mode;
end