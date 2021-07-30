%---------------evelope function by hilbert function----------
function [env,pha]=envelope(dat0)
    hilb=hilbert(dat0);
    env=abs(hilb);
    pha=atan2(imag(hilb),real(hilb));
    return
    %---------------copy from Chuck's "surface.m" program---------
    hilb=hilbert(dat0);
    ihilb=imag(hilb);
    env=abs(dat0 - i*ihilb);
    pha=atan2(-ihilb,dat0);
    return
