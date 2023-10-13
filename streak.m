

clc
clear



add2='FLASH.DIAG/CAMERA/OTR9FL2XTDS/SPECTRUM.X.SIG';
%    'FLASH.DIAG/CAMERA/OTR9FL2XTDS/SPECTRUM.X.CM'




for j=1:10
    x(1,j)=getfield(doocsread(add2), 'data');
    pause(0.12)

end

add_phase='FLASH.RF/LLRF.CONTROLLER/CTRL.POLARIX/SP.PHASE';
phase=getfield(doocsread(add_phase), 'data');


tem=doocswrite(add_phase, phase-180)

pause(1)



for j=1:10
    x(3,j)=getfield(doocsread(add2), 'data');
    pause(0.12)

end



addr_xtds_onoff='FLASH.DIAG/TIMINGINFO/FLFXTDS/ON_BEAM'
tmp      = doocswrite(addr_xtds_onoff, 0);

pause(1)
for j=1:10
    x(2,j)=getfield(doocsread(add2), 'data');
    pause(0.12)

end


tmp      = doocswrite(addr_xtds_onoff, 1);
tem=doocswrite(add_phase, phase)

%%
xav=mean(x,2)
%%
figure; plot([-1,0,1], xav, '-o')