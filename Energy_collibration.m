
clear all 
add1='FLASH.RF/LLRF.SUMVOLTAGE_CTRL/ACC45/BR2.SP.SUMVOLTAGE';
t=doocsread(add1)
volt=t.data

%%
%add2='FLASH.DIAG/CAMERA/OTR9FL2XTDS/SPECTRUM.Y.TD'
%add2='FLASH.DIAG/CAMERA/OTR9FL2XTDS/SPECTRUM.Y.MEAN';
add2='FLASH.DIAG/CAMERA/OTR9FL2XTDS/SPECTRUM.Y.CM'
%y=doocsread(add2);
%y.data

%%
i=0;
for v=volt-2:1: volt+2
    i=i+1;
    tem=doocswrite(add1, v)
    for j=1:20
    x(i,j)=getfield(doocsread(add1), 'data');
    Center=doocsread(add2);
    y(i,j)=Center.data;
    pause(0.1)
    end
    pause(0.2)

end
tem=doocswrite(add1, volt)
%%

figure ()
plot(x,y, 'o')
hold on 
plot(mean(x ,2), mean(y , 2))

%%
name='energy_coll'
energy_coll=struct
energy_coll.volt=x;
energy_coll.y=y;

save('energy_coll', 'energy_coll', '-mat')
