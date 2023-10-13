clear all
x=10
% feedback
add_COFM='FLASH.DIAG/CAMERA/OTR9FL2XTDS/SPECTRUM.X.CM'
add2='FLASH.DIAG/CAMERA/OTR9FL2XTDS/SPECTRUM.X.MEAN';
add_phase='FLASH.RF/LLRF.CONTROLLER/CTRL.POLARIX/SP.PHASE'
while x>1
    for j=1:20
        Center=doocsread(add_COFM);
        x(j)=Center.data;
        pause(0.12)
    end
    phase_XTDS=getfield(doocsread(add_phase), 'data');

    if mean(x) > 8
        tem=doocswrite(add_phase, phase_XTDS+1)
    end
    if mean(x) < 5
        tem=doocswrite(add_phase, phase_XTDS-1)
    end

end
